#using AmbulanceDeployment
using DataFrames, Winston, JLD, CSV, Gurobi, JuMP

#HOW DATAFRAMES WORK
#https://julia.school/julia/dataframes/#what-are-dataframes

hourly_calls = CSV.File("..\\test\\Austin_Data\\Full_WeekdayCalls.csv") |> DataFrame
adjacent_nbhd = CSV.File("..\\test\\Austin_Data\\adjacent_nbhd.csv") |> DataFrame
coverage = CSV.read("..\\test\\Austin_Data\\coverage_real.csv", DataFrame, header=false)
incidents = CSV.File("..\\test\\Austin_Data\\austin_incidents.csv") |> DataFrame

regions = Int[parse(Int,string(x)) for x in names(hourly_calls[:,6:ncol(hourly_calls)])]
locations = collect(1:size(coverage,2))
adjacent = convert(Array, adjacent_nbhd[:,2:ncol(adjacent_nbhd)])[regions,regions] .> 0.5
demand = convert(Array,hourly_calls[:,6:end]);

regions2index = Dict{Int,Int}(regions[i]=>i for i in 1:length(regions))
#= previously there was an error that didnt account for all the regions because it
 it only included regions for incidents. the regions without incidents are changed
 to 0=#
for x in incidents[!,:neighborhood]
    if(!haskey(regions2index,x))
        regions2index[x] = 0
    end
end

# We focus on emergency calls during the "peak period" (8AM - 8PM),
# with the emergency calls from the first 3 month as our training set,
# and the subsequent emergency calls from the remaining months as our test set
peak_period = (hourly_calls[!,:hour] .>= 8) .* (hourly_calls[!,:hour] .<= 20) #indicator vector #half of 10k calls are in peak period
indices = 1:DataFrames.nrow(hourly_calls);
train_filter = (hourly_calls[!,:year] .== 2019) .* (hourly_calls[!,:month] .<= 3) #indicator function... trash one
test_filter  = .~train_filter;

train_indices = indices[train_filter]
test_indices = indices[test_filter]
debug = true;
if(debug == true)
    debug_length = 20
    train_indices = train_indices[1:debug_length]
    test_indices = test_indices[1:debug_length]
end

namb = 37
p = DeploymentProblem(namb, length(locations), length(regions), demand, train_indices,
test_indices, coverage[regions,:], Array{Bool,2}(adjacent));


amb_deployment = Dict{Symbol, Dict{Int, Vector{Int}}}()
#amb_deployment = Dict{String, deployment_model}()

# for (next_deployment_model, name) in ((next_dp -> StochasticDeployment(next_dp, nperiods=500), :Stochastic),
#                                 (next_dp -> MEXCLPDeployment(next_dp, 0.654), :MEXCLP),
#                                 (next_dp -> MALPDeployment(next_dp, 0.654), :MALP))
#    println("$name: ")
#    amb_deployment[name] = Dict{Int, Vector{Int}}()
#    for namb in 30:5:50
#        println("$namb ")
#        p.nambulances = namb
#        next_model = next_deployment_model(p)
#        set_optimizer(next_model.m, Gurobi.Optimizer)
#        @time optimize!(next_model, p)
#        amb_deployment[name][namb] = deployment(next_model)
#    end
#    println()
# end

#This creates a symbol, which seems to be a fancy string.
# #https://stackoverflow.com/questions/23480722/what-is-a-symbol-in-julia
# (deployment_model, name) = (next_dp -> StochasticDeployment(next_dp, nperiods=500), :Stochastic)
#
# model = StochasticDeployment(p, nperiods=500)
# set_optimizer(model.m, Gurobi.Optimizer)
# @time optimize!(model, p)
# amb_deployment[name][namb] = deployment(model)

JLD.jldopen("team_stats.jld", "w") do file
    write(file, "amb_deployment", amb_deployment)
end

#check feasibility
I = 1:p.nlocations
J = 1:p.nregions
for j in J
    count = 0;
    for i in filter(i->p.coverage[j,i], I)
        count = count + 1;
    end
    if(count == 0)
        print(count)
    end
end
