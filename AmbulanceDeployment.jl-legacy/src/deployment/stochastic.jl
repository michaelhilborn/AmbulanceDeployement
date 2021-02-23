#=
Author : Ng Yeesian
Modified : Guy Farmer
generates the stochastic deployment model
=#
struct StochasticDeployment <: DeploymentModel
    m::JuMP.Model
    x::Array{JuMP.VariableRef,1}
    y::Array{JuMP.VariableRef,3}
end
deployment(m::StochasticDeployment) = [round(Int,x) for x in JuMP.value.(m.x)]

function StochasticDeployment(p::DeploymentProblem; nperiods=params.nperiods, tol=params.δ,
    #solver=GurobiSolver(OutputFlag=0))
    solver=Gurobi.Optimizer(OutputFlag=0))

    nperiods = min(length(p.train), nperiods)
    demand = p.demand[p.train,:]
    I = 1:p.nlocations
    J = 1:p.nregions
    T = 1:nperiods

    m = Model(GLPK.Optimizer)
    JuMP.@variable(m, x[1:p.nlocations] >= 0, Int)
    JuMP.@variable(m, y[1:p.nlocations,1:p.nregions,1:nperiods] >= 0, Int)
    JuMP.@variable(m, z[1:p.nregions,1:nperiods] >= 0, Int)

    JuMP.@objective(m, Min, sum(z[j,t] for j=J, t=T) + tol*sum(y[i,j,t] for i=I, j=J, t=T))

    JuMP.@constraint(m, sum(x[i] for i=I) <= p.nambulances)

    for j in J # coverage over all regions
        JuMP.@constraint(m, sum(x[i] for i in filter(i->p.coverage[j,i], I)) >= 1)
    end

    # flow constraints at each station
    for i in I, t in T
        outflow = JuMP.@expression(m, sum(y[i,j,t] for j in filter(j->p.coverage[j,i], J)))
        JuMP.@constraint(m, x[i] >= outflow)
    end

    # shortfall from satisfying demand/calls
    for j in J, t in T
        inflow = JuMP.@expression(m, sum(y[i,j,t] for i in filter(i->p.coverage[j,i], I)))
        JuMP.@constraint(m, z[j,t] >= demand[t,j] - inflow)
    end

    for i in I
        JuMP.@constraint(m, x[i] <= 5)
    end

    StochasticDeployment(m, x, y) #y
end

optimize!(model::StochasticDeployment) = JuMP.optimize!(model.m)
optimize!(model::StochasticDeployment, p::DeploymentProblem) = JuMP.optimize!(model.m)
