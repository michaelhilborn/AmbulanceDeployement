#=
Author : Ng Yeesian
Modified : Guy Farmer / Michael Hilborn
generates the assignment model / assigns ambulances based on a queueing system
currently operates on gurobi 0.8.1
=#
struct AssignmentModel <: RedeployModel
    # model::Gurobi.Model
    model::JuMP.Model
    lambda::Float64

    hosp2stn::Matrix{Float64}
    stn2stn::Matrix{Float64}

    assignment::Vector{Int} # which location the ambulance is assigned to
    ambulances::Vector{Vector{Int}} # list of ambulances assigned to each location
    status::Vector{Symbol} # the current status of the ambulance
        # possible statuses: :available, :responding, :atscene, :conveying, :returning
    fromtime::Vector{Int} # the time it started the new status
    hospital::Vector{Int} # the hospital the ambulance is at (0 otherwise)

    soln::Vector{Float64} # buffer for storing dynamic assignment solutions
end

function AssignmentModel(
        p::DeploymentProblem,
        available::Vector{Int},
        # utilization::Vector{Float64},
        hospitals::DataFrame,
        stations::DataFrame;
        lambda::Float64 = 100.0
    )

    nambulances = sum(available)
    nlocations = length(available)
    assignment = zeros(Int, nambulances)

    hosp2stn = convert(Matrix{Float64},
        hcat([hospitals[Symbol("stn$(i)_min")] for i in 1:nlocations]...)
    )
    stn2stn =  convert(Matrix{Float64},
        hcat([stations[Symbol("stn$(i)_min")] for i in 1:nlocations]...)
    )

    k = 1
    ambulances = [Int[] for i in 1:nlocations]
    for i in eachindex(available), j in 1:available[i]
        assignment[k] = i
        push!(ambulances[i], k)
        k += 1
    end
    @assert k == nambulances + 1

    @assert sum(length(a) for a in ambulances) == nambulances
    status = fill(:available, nambulances)
    fromtime = zeros(Int, nambulances)
    hospital = zeros(Int, nambulances)

    #m = Gurobi.Model(Gurobi.Env(), "redeploy", :minimize)
    m = Model(GLPK.Optimizer)


    # inquire what is the purpose of the w variable for simulation
    # guess - since its a binary variable of ambulances and locations, a 1 means youre sending a specific amb to a specific location
    # its a weight. weight 1 of amb 1 plus weight 2 of amb 2 plus weight 3 of amb 3 = 1
    # inquire what specifically is lambda, and what is the difference between 0 and 100
    #Gurobi.setparam!(m, "OutputFlag", 0)
    @variable(m, w[1:nambulances, 1:nlocations], Bin) # w variable

    @variable(m, eta1[1:nlocations], lower_bound = -Inf, upper_bound = Inf)# eta1
    @variable(m, eta2[1:nambulances, 1:nlocations], lower_bound = 0, upper_bound = Inf) # eta2
    @variable(m, eta3[1:nlocations], lower_bound = 0., upper_bound = Inf)# eta3 := eta1^2


    # η₁[i] >= available[i] - sum(w[a,i] for a in 1:nambulances)
    #     reformulated to
    # sum(w[a,i] for a in 1:nambulances) + η₁[i] >= available[i]
    # for i in 1:nlocations
    #         Gurobi.add_constr!(m,[((1:nambulances).-1)*nlocations.+i; nambulances*nlocations + i], ones(nambulances + 1), '>', Float64(available[i]))
    #        end
    for i in 1:nlocations
        @constraint(m, sum(w[a,i] for a in 1:nambulances) + 1 >= x[i])
    end

    ## repaired with dot syntax

    # sum(w[a,i] for i in 1:nlocations) == 1       [a=1:nambulances]
    # for a in 1:nambulances
    #     Gurobi.add_constr!(m,
    #         collect((a-1)*nlocations .+ (1:nlocations)), # inds
    #         ones(nlocations), # coeffs
    #         '=', 1.)
    # end

    for a in 1:nambulances
           @constraint(m, sum(w[a,i] for i in nlocations) == 1)
       end

    # eta2[a,i] >= |w[a,i] - (assignment[a] == i)|   [a=1:nambulances, i=1:nlocations]
    #     reformulated to
    # eta2[a,i] >= w[a,i] - (assignment[a] == i)   i.e. eta2[a,i] - w[a,i] >= -(assignment[a] == i)
    # eta2[a,i] >= - w[a,i] + (assignment[a] == i) i.e. eta2[a,i] + w[a,i] >=  (assignment[a] == i)
    for a in 1:nambulances, i in 1:nlocations
        offset = (a-1)*nlocations + i
        inds = [(nambulances+1)*nlocations + offset, offset]
        # add_constr(model, left hand side of equation, separator, right hand side of equation)
        Gurobi.add_constr!(m, inds, [1., -1.], '>', - Float64(assignment[a] == i))
        Gurobi.add_constr!(m, inds, [1., 1.], '>', Float64(assignment[a] == i))
    end

    # eta3[i] = eta1[i]^2
    #     reformulated into
    # eta3[i] >= 0
    # eta3[i] >= 1 + 2*eta1[i]
    # eta3[i] >= 4 + 4*eta1[i]
    # eta3[i] >= 9 + 6*eta1[i]
    # ...
    for i in 1:nlocations
        for k in 0.1:0.1:0.9
            Gurobi.add_constr!(m,
                [2*nambulances*nlocations + nlocations + i, nambulances*nlocations + i], # inds
                [1., -2*k], # coeffs
                '>', Float64(k^2))
        end
        for k in 1:3
            Gurobi.add_constr!(m,
                [2*nambulances*nlocations + nlocations + i, nambulances*nlocations + i], # inds
                [1., -Float64(2*k)], # coeffs
                '>', Float64(k^2))
        end
    end

    AssignmentModel(m, lambda, hosp2stn, stn2stn, assignment, ambulances,
                    status, fromtime, hospital, zeros(nambulances*nlocations))
end

function reassign_ambulances!(
        ems,
        problem::DispatchProblem,
        redeploy::AssignmentModel,
        t::Int
    )
    nlocations = length(redeploy.ambulances)
    nambulances = length(redeploy.assignment)
    # DEBUG
    # for a in 1:nambulances
    #     @assert t >= redeploy.fromtime[a]
    # end
    cost(a,i) = if redeploy.status[a] == :available
             0                                 + redeploy.stn2stn[redeploy.assignment[a],i]
        elseif redeploy.status[a] == :responding
            55 - (t - redeploy.fromtime[a])/60 + redeploy.stn2stn[redeploy.assignment[a],i] # + redeploy.hosp2stn[redeploy.hospital[a],i]
        elseif redeploy.status[a] == :atscene
            45 - (t - redeploy.fromtime[a])/60 + redeploy.stn2stn[redeploy.assignment[a],i] # + redeploy.hosp2stn[redeploy.hospital[a],i]
        elseif redeploy.status[a] == :conveying
            30 - (t - redeploy.fromtime[a])/60 + redeploy.stn2stn[redeploy.assignment[a],i] # + redeploy.hosp2stn[redeploy.hospital[a],i]
        elseif redeploy.status[a] == :returning
            15 - (t - redeploy.fromtime[a])/60 + redeploy.stn2stn[redeploy.assignment[a],i]
        elseif redeploy.status[a] == :redeploying
            10 - (t - redeploy.fromtime[a])/60 + redeploy.stn2stn[redeploy.assignment[a],i]
        end
    # (1) Optimize Dynamic Assignment Problem
    # (1i) change RHS to account for wait_queue
    nwait = [length(problem.wait_queue[nbhd]) for nbhd in 1:size(problem.coverage,1)]
    Gurobi.set_dblattrarray!(redeploy.model, "RHS", 1, nlocations, [
        Float64(problem.deployment[i] + sum(
            nwait[nbhd] for nbhd in 1:size(problem.coverage,1)
            if problem.coverage[nbhd,i]
        ))
        for i in 1:nlocations
    ])
    # (1ii) change coeffs to account for traveling time
    let con = Cint[0], ind = Cint[0], val = Float64[0.0]
        for i in 1:nlocations
            con[1] = i
            for a in 1:nambulances
                ind[1] = (a-1)*nlocations + i
                val[1] = 1-max(min(60,cost(a,i)),0)/45
                Gurobi.chg_coeffs!(redeploy.model, con, ind, val)
            end
        end
    end
    Gurobi.update_model!(redeploy.model)
    Gurobi.optimize(redeploy.model)
    # @show Gurobi.get_objval(redeploy.model)

    # (2) Reassign ambulances based on optimal solution
    Gurobi.get_dblattrarray!(redeploy.soln, redeploy.model, "X", 1)
    for a in 1:nambulances
        let stn = redeploy.assignment[a]
            for i in 1:nlocations
                if stn != i && redeploy.soln[(a-1)*nlocations + i] > 0.5
                    # redeploy an existing ambulance
                    push!(problem.redeploy_events, (a,redeploy.assignment[a],i,t))
                    if redeploy.status[a] == :available
                        @assert problem.available[redeploy.assignment[a]] > 0
                        problem.available[stn] -= 1
                        t_end = t + ceil(Int, 0*60*redeploy.stn2stn[stn,i])
                        redeploying_to!(redeploy, a, i, t)
                        enqueue!(ems.eventqueue, (:done, 0, t_end, a), t_end)
                    else
                        println("redeploying amb $a from $(redeploy.assignment[a]) to $i")
                        ambulances = redeploy.ambulances[redeploy.assignment[a]]
                        @assert !in(a, redeploy.ambulances[i])
                        @assert !in(a, ambulances)
                        in(a, ambulances) && deleteat!(ambulances, findfirst(ambulances, a))
                        @assert !in(a, ambulances)
                        redeploy.assignment[a] = i
                    end
                end
            end
        end
    end
end
