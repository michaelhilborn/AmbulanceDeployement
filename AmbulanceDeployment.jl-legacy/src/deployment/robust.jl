
include("../problem.jl")
import JuMP, JuMP.optimize!
Pkg.add("JuMP")
Pkg.add("GLPK")
#Pkg.add(name=”Gurobi” version=”0.8.1”)
using JuMP, Gurobi, GLPK, Distributions

struct Gamma
    _single::Vector{Int}
    _local::Vector{Int}
    _regional::Vector{Int}
    _global::Int
end
#remember to comment out with package
struct Params
    α::Float64 # Probabilistic Guarantee
    ε::Float64 # Convergence
    δ::Float64 # Solver Tolerance

    nperiods::Int # for StochasticDeployment

    maxiter::Int # for RobustDeployment
end
paramss = Params(0.01, 0.5, 1e-6, 500, 50)
#might need to change JuMP.Variable into JuMP.VariableRef
struct Qrobust
    m::JuMP.Model
    I::UnitRange{Int}
    J::UnitRange{Int}
    d::Vector{JuMP.VariableRef}
    q::Vector{JuMP.VariableRef}
    γ::Gamma
end

struct RobustDeployment <: DeploymentModel
    m::JuMP.Model
    Q::Qrobust

    I::UnitRange{Int}
    J::UnitRange{Int}

    x::Vector{JuMP.VariableRef}
    y::Vector{Matrix{JuMP.VariableRef}}
    z::Vector{Vector{JuMP.VariableRef}}
    η::JuMP.VariableRef

    scenarios::Vector{Vector{Int}}
    upperbounds::Vector{Float64}
    lowerbounds::Vector{Float64}
    deployment::Vector{Vector{Int}}

    upptiming::Vector{Float64}
    lowtiming::Vector{Float64}
end
deployment(m::RobustDeployment) = m.deployment[end]

function Gamma(p::DeploymentProblem; α=paramss.α)
    demand = p.demand[p.train,:]
    γ_single = Any[] # vec(maximum(demand,1) + 1 *(maximum(demand,1) .== 0))

    for x in 1:size(demand)[2]
        y = maximum(demand[:,x])
        if y == 0
            y = 1
        end
        push!(γ_single,y)
    end
#have to know what the purpose is for dealing with the mean in this manner,
# quick repair is to adjust the dimensions by changing sum(demand[:,vec(p.adjacency[i,:])], 2) to sum(demand[:,vec(p.adjacency[i,:])])
# this will give the total sum across 2 dimensions, but will nullify the mean function
# will change the mean function from mean(A) to actually divide the sum by the total indices in the demand subarray
# the same logic is occuring for γ_regional and γ_global
# /90 represents the number of calls in a 3 month period
    #γ_local = [quantile(Poisson((sum(demand[:,vec(p.adjacency[i,:])], 2))/90),1-α) for i=1:p.nregions]
    # γ_local = [quantile(Poisson(mean(sum(demand[:,vec(p.adjacency[i,:])]))),1-α) for i=1:p.nregions]
    γ_local = [quantile(Poisson(mean(sum(demand[:,vec(p.adjacency[i,:])]))/90),1-α) for i=1:p.nregions]
    γ_regional = [quantile(Poisson(mean(sum(demand[:,p.coverage[:,i]]))/90),1-α) for i in 1:p.nlocations]
    γ_global = quantile(Poisson(mean(sum(demand))/90),1-α)
    Gamma(γ_single,γ_local,γ_regional,γ_global)
end

function Qrobust(problem::DeploymentProblem; α=paramss.α, verbose=false, solver=Gurobi.Optimizer(OutputFlag=0)) #, MIPGapAbs=0.9)
    if verbose
        solver=Gurobi.Optimizer(OutputFlag=1) #, MIPGapAbs=0.9)
        #solver=JuMP.set_optimizer_attribute(model.m,MOI.TimeLimitSec(), 30)
    end
    γ = Gamma(problem, α=α)
    upp_bound = maximum(γ._single)
    I = 1:problem.nlocations
    J = 1:problem.nregions

    m = Model(GLPK.Optimizer)
    JuMP.@variable(m, d[1:problem.nregions]>=0, Int)
    JuMP.@variable(m, p[1:problem.nregions], Bin)
    JuMP.@variable(m, q[1:problem.nlocations], Bin)

    for i in I, j in J
        problem.coverage[j, i] && JuMP.@constraint(m, p[j] <= q[i])
    end

    # Uncertainty
    for j in J
        JuMP.@constraint(m, d[j] <= γ._single[j]*p[j])
        adjacent_regions = filter(k->problem.adjacency[k,j],J)
        JuMP.@constraint(m, sum(d[k] for k in adjacent_regions) <= γ._local[j])
    end
    for i in I
        covered_regions = filter(j->problem.coverage[j,i],J)
        JuMP.@constraint(m, sum(d[j] for j in covered_regions) <= γ._regional[i])
    end
    JuMP.@constraint(m, sum(d[j] for j in J) <= γ._global)

    Qrobust(m, I, J, d, q, γ)
end

function evaluate(Q::Qrobust, x::Vector{T}) where {T <: Real}
    JuMP.@objective(Q.m, Max, sum(Q.d[j] for j in Q.J) - sum(x[i]*Q.q[i] for i in Q.I))
    #status = JuMP.solve(Q.m)
    status = JuMP.optimize!(Q.m)
    JuMP.objective_value(Q.m), Int[round(Int,JuMP.value(d)) for d in Q.d]

end

function evaluate_objvalue(Q::Qrobust, x::Vector{T}) where {T <: Real}
    JuMP.@objective(Q.m, Max, sum(Q.d[j] for j in Q.J) - sum(x[i]*Q.q[i] for i in Q.I))
    status = JuMP.solve(Q.m)
    JuMP.objective_value(Q.m)
end


function RobustDeployment(p::DeploymentProblem; α=paramss.α)
     eps=paramss.ε
     tol=paramss.δ
     m = Model(GLPK.Optimizer)
     solver=Gurobi.Optimizer(OutputFlag=0, MIPGapAbs=0.9)
     #JuMP.set_optimizer_attribute(m,MOI.TimeLimitSec(), 30)
     verbose=false
     master_verbose=false
    if master_verbose
        solver=Gurobi.Optimizer(OutputFlag=1, MIPGapAbs=0.9)
        #JuMP.set_optimizer_attribute(m,MOI.TimeLimitSec(), 30)
    end
    I = 1:p.nlocations
    J = 1:p.nregions

    warmstart = naive_solution(p)

    #m = JuMP.Model(solver=solver)
    # m = Model(with_optimizer(GLPK.Optimizer))
    JuMP.@variable(m, x[i=1:p.nlocations] >= 0, Int, start=warmstart[i])
    JuMP.@variable(m, η >= 0)
    y = Vector{Matrix{JuMP.VariableRef}}()
    z = Vector{JuMP.VariableRef}()

    # Initial Restricted Master Problem
    JuMP.@objective(m, Min, η)
    JuMP.@constraint(m, sum(x[i] for i=I) <= p.nambulances)
    for j in J # coverage over all regions
        JuMP.@constraint(m, sum(x[i] for i in filter(i->p.coverage[j,i], I)) >= 1)
    end

    RobustDeployment(m, Qrobust(p, α=α, verbose=verbose), I, J, x, y, z, η,
                     Vector{Vector{Int}}(), Vector{Float64}(), Vector{Float64}(),
                     Vector{Int}[warmstart], Vector{Float64}(), Vector{Float64}())
end

function add_scenario(model::RobustDeployment, p::DeploymentProblem, scenario::Vector{T}; tol=paramss.δ) where {T <: Real}
    # Create variables yˡ
    push!(model.y, Array{VariableRef}(undef,(p.nlocations,p.nregions)))
    l = length(model.y)
    # for i in model.I, j in model.J
    #     model.y[l][i,j] = JuMP.VariableRef(model.m, lower_bound=0, upper_bound=p.nambulances, variable_type=:Int, base_name=String("y[$i,$j,$l]"))
    # end
    for i in model.I, j in model.J
               variable = @variable(model.m)
               set_lower_bound(variable, 0)
               set_upper_bound(variable, p.nambulances)
               set_integer(variable)
               set_name(variable, String("y[$i,$j,$l]"))
               model.y[l][i,j] = variable
           end
    push!(model.z, Array{VariableRef}(undef, p.nregions))
    # for j in model.J
    #     model.z[l][j] = JuMP.VariableRef(model.m, 0, Inf, :Int, String("z[$j,$l]"))
    # end
    for j in model.J
               variable = @variable(model.m)
               set_lower_bound(variable, 0)
               set_upper_bound(variable, Inf)
               set_integer(variable)
               set_name(variable, String("z[$j,$l]"))
               model.z[l][j] = variable
           end

    # (1) η >= 1ᵀ(dˡ + Bᴶyˡ)^+
    JuMP.@constraint(model.m, model.η >= sum(model.z[l][j] for j=model.J) + tol*sum(model.y[l][i,j] for i=model.I, j=model.J))
    for i in model.I
             outflow = JuMP.@expression(model.m, sum(model.y[l][i,j] for j in filter(j->p.coverage[j,i], model.J)))
             JuMP.@constraint(model.m, model.x[i] >= outflow)
           end
    # (2) yˡ ∈ Y(x)
    for j in model.J # shortfall from satisfying demand/calls
        inflow = JuMP.@expression(model.m, sum(model.y[l][i,j] for i in filter(i->p.coverage[j,i], model.I)))
        JuMP.@constraint(model.m, model.z[l][j] >= scenario[j] - inflow)
    end
end

function optimize!(model::RobustDeployment, p::DeploymentProblem; verbose=false, maxiter=params.maxiter, eps=params.ε)
    LB = 0.0
    UB, scenario = evaluate(model.Q, model.deployment[end])
    push!(model.lowerbounds, LB)
    push!(model.upperbounds, UB)
    push!(model.scenarios, scenario)

    for k in 1:maxiter
        verbose && println("iteration $k: LB $LB, UB $UB")
        abs(UB - LB) < eps && break
        verbose && println("  solving Q with $(model.deployment[end])")

        add_scenario(model, p, scenario)
        #tic()
        #status = JuMP.solve(model.m)
        JuMP.optimize!(model.m)
        #push!(model.lowtiming, toq())
        #@assert status == :Optimal

        LB = JuMP.objective_value(model.m)

        push!(model.deployment, [round(Int,x) for x in JuMP.value.(model.x)])

        #tic()
        shortfall, scenario = evaluate(model.Q, model.deployment[end])
        #push!(model.upptiming, toq())
        UB = min(UB, shortfall)

        # for tracking convergence later
        push!(model.upperbounds, UB)
        push!(model.scenarios, scenario)
        push!(model.lowerbounds, LB)
    end
end

optimize!(model::RobustDeployment, p::DeploymentProblem) = JuMP.optimize!(model.m)
