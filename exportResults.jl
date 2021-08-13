"""
export results of the model

beware that some field may be changed according to the experiment

"""


# export lines in a .csv file
df3 = DataFrame(Source = Int[], Sink = Int[], Route = Int[], flow = Float64[])
for i in N, j in N, t in T
    if pth[i,j] > 0 && value(d[i, j, t]) > 0.0000001
        push!(df3, (i, j, t, value(d[i,j,t])))
    end
end
CSV.write("/Users/akrm/Dropbox/Julia/transitComm/encapsulated/output/routes.csv",  DataFrame(df3), writeheader=true)

# store the results from problem
# the objective value
pTranCurrent = objective_value(pTran)

# export passenger flows in a .csv file
df = DataFrame(Source = Int[], Sink = Int[], Commodity = Int[], Route = Int[], flow = Float64[], distance = Int64[])
for i in N, j in N, k in N, t in T
    if pth[i,j] > 0 && value(x[i, j, k, t]) > 0
        push!(df, (i, j, k, t, value(x[i,j,k,t]), pth[i,j]))
    end
end
CSV.write("/Users/akrm/Dropbox/Julia/transitComm/encapsulated/output/solutions.csv",  DataFrame(df), writeheader=true)

# export frequencies in a .csv file
df2 = DataFrame(freq = Float64[])
for t in T
        push!(df2, value(f[t]))
end
CSV.write("/Users/akrm/Dropbox/Julia/transitComm/encapsulated/output/freqs.csv", DataFrame(df2), writeheader = true)


df5 = DataFrame(Source=Int[], Sink=Int[], Commodity=Int[], Line=Int[], Flow=Float64[], distance = Int64[])
for i in N, j in destination[i], k in passenger[i], t in T
    if value(trnsfer[i,j,k, t]) > 0
        push!(df5, (i,j,k,t,value(trnsfer[i,j,k,t]), pth[i,j]))
    end
end
CSV.write("/Users/akrm/Dropbox/Julia/transitComm/encapsulated/output/penalty.csv", DataFrame(df5), writeheader=true)

#
df4 = DataFrame(Source = Int[], Route = Int[], val = Float64[])
for i in N, t in T
    if value(u[i, t]) > 0
        push!(df4, (i, t, value(u[i,t])))
    end
end
CSV.write("/Users/akrm/Dropbox/Julia/transitComm/encapsulated/output/orderoftheline.csv",  DataFrame(df4), writeheader=true)

#=
df5= DataFrame(Source = Int[], Sink = Int[], Line = Int[], Value = Float64[])
for i in N, j in N, t in T
    if pth[i,j]>0 && value(y[i, j, t]) >=0
        push!(df5, (i, j, t, value(y[i,j,t])))
    end
end
CSV.write("/Users/akrm/Dropbox/Julia/transitComm/encapsulated/output/comSubtour.csv",  DataFrame(df4), writeheader=true)
=#

#= export the routes for d variable with 2 index
df3 = DataFrame(Source = Int[], Route = Int[], flow = Float64[])
for i in N, t in T
    if value(d[i, t]) > 0.001
        push!(df3, (i, t, value(d[i,t])))
    end
end
CSV.write("/Users/akrm/Dropbox/Julia/transitComm/encapsulated/output/routesD-2INDEX.csv",  DataFrame(df3), writeheader=true)
=#
