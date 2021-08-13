"""
the mathematical model with MTZ constraints
the model is the rigorous one
no experiment, it is just the model
"""

flowq = ODistance(true,false)
pth = ODistance(false,true)
sup, rhs = rightHandSide()
params = problemParameters()
N, NT, T, destination, passenger,lineDict = requiredDictionaries()

# mathematical model
pTran = Model(optimizer_with_attributes(Gurobi.Optimizer,
                "Presolve" => 0,
                "MIPGap" => 0.0001, "MIPFocus" => 1,
                #"Heuristics" => 0.95, "Cuts" => 0, "ImproveStartTime" => 400,
                #"Threads" => 4,
                #"Method" => 1,
                #"MarkowitzTol" => 0.2,
                #"VarBranch" => 3,
                #"RINS" => 100,
                #"IntFeasTol" => 10e-3,
                #"VarBranch" => 2,
                "TimeLimit" => 1000
                ))

# add norelheur heuristics in which heusristics are used before relaxation for the cases relaxation takes time
#set_optimizer_attribute(pTran,"NoRelHeurTime",500)

# passenger flow decision variable
@variable(pTran, x[i in N, j in destination[i], k in passenger[i], t in T], lower_bound=0)

# transfer flow decision variable
@variable(pTran, trnsfer[i in N, j in destination[i], k in passenger[i], t in T], lower_bound=0)

# line flow decision variable
@variable(pTran, d[i in NT, j in lineDict[i], t in T], Bin)

# adding variable frequency of lines
@variable(pTran, f[t in T], lower_bound=0.0)

# adding MTZ subtour elimination variable
@variable(pTran, u[i in N, t in T], lower_bound=0.0)

# LINE FLOW CONSERVATION CONSTRAINTS:
# lines super source node sends 1 unit for each line
@constraint(pTran, lineBal[t in T],
                            sum(d[params.nNodes+1,j,t] for j in lineDict[params.nNodes+1]) == 1)

# line flows pass through stop nodes
@constraint(pTran, lineBal2[i in N, t in T],
                            + sum(d[i,j,t] for j in lineDict[i])
                            - sum(d[j,i,t] for j in destination[i])
                            - d[params.nNodes+1, i, t] == 0)

# line super sink node's flow demand is satisfied
@constraint(pTran, lineBal3[t in T],
                            - sum(d[i,params.nNodes+2,t] for i in N) == -1)
#

#terminals = [2,7,9,12]
# only potential terminal set can be a terminal
# @constraint(pTran, lineBal4[i in N,t in T; i ∉ terminals], d[params.nNodes+1, i, t] <= 0)
# @constraint(pTran, lineBal5[i in N,t in T; i ∉ terminals], d[i, params.nNodes+2, t] <= 0)

# line directions are already works bidirectional
@constraint(pTran, bidirect[i in N, j in destination[i], t in T],
                            d[i,j,t] + d[j,i,t] <= 1)

# only one unit of line flow outs from a node in a line
@constraint(pTran, outFlow[i in N, t in T],
                            sum(d[i,j,t] for j in destination[i]) <= 1)

# only one unit of line flow outs from a node in a line
@constraint(pTran, outFlow2[j in N, t in T],
                            sum(d[i,j,t] for i in destination[j]) <= 1)
#=
@constraint(pTran, traffic[i in N, j in destination[i]],
                            sum(d[i,j,t] for t in T) <= 2)
=#

# PASSENGER FLOW CONSERVATION CONSTRAINTS:
# supply of passenger flows
@constraint(pTran, flowb[i in N],
                            sum(x[i,j,i,t] for j in destination[i], t in T) == rhs[i,i])

# passenger flows leaves lines at the desired stop, otherwise just pass over that stops
@constraint(pTran, flowb2[i in N, k in passenger[i]; k!=i],
                            + sum(x[i,j,k,t] for j in destination[i], t in T)
                            - sum(x[j,i,k,t] for j in destination[i], t in T) == rhs[i,k])

# passenger arrives a stop from the previous stop of the line or
# transfer from another line (except the stop's own passengers)
@constraint(pTran, transfer[i in N, j in destination[i], k in passenger[i], t in T; (k!=i)],
                            x[i,j,k,t] <= sum(x[r,i,k,t] for r in destination[i] if r!=j)
                                                + trnsfer[i,j,k,t])
#
@constraint(pTran, transfer2, sum(trnsfer[i,j,k,t] for i in N, j in destination[i], k in passenger[i], t in T)
                                <= sum(sup)*0.083)
#

#=
@constraint(pTran, transfer3, sum(trnsfer[i,j,k,t] for i in N, j in destination[i], k in passenger[i], t in T)
                                >= sum(sup)*0.08)
 =#

# a passenger cannot return its own supply stop
# @constraint(pTran, ownCom[i in N, j in destination[i], t in T], x[i,j,j,t] <= 0)

# a passenger can through one stop to another iff there is line flow over there
# in other words the segment of a line is defined
@constraint(pTran, capacity[i in N, j in destination[i], t in T],
                            sum(x[i,j,k,t] for k in passenger[i])
                                <=params.bigNumber*params.cap[t]*(d[i,j,t] + d[j,i,t]))

#@constraint(pTran, capacity2[i in N, j in destination[i], k in passenger[i], t in T; i<j],
#                            x[i,j,k,t] + x[j,i,k,t] <= sup[k]*d[i,j,t] + d[j,i,t])


# frequency determination constraint
@constraint(pTran, freq[i in N, j in destination[i], t in T],
                           sum(x[i,j,k,t] for k in passenger[i]) + sum(x[j,i,k,t] for k in passenger[j])
                                <= f[t]*params.cap[t])

# limits the number of stops in a line
@constraint(pTran, limitItinerary[t in T],
                            sum(d[i,j,t] for i in N, j in destination[i]) == params.itinerary-1)
#

# subtour elimination by MTZ
@constraint(pTran, tour[i in N, j in destination[i], t in T],
                            u[i,t] - u[j,t] + (params.itinerary-1)*(d[i,j,t]) <= params.itinerary-2)

@constraint(pTran, tour2[i in N, t in T],
                            u[i,t] <= params.itinerary-1)
#

# a stop must be located in at least one line
#@constraint(pTran, locat[i in N], sum(d[i,j,t] for j in destination[i], t in T) >= 1)

# there will be a line segment iff there is a passenger flow
# @constraint(pTran, locat2[i in N, j in destination[i], t in T], d[i,j,t] <= sum(x[i,j,k,t] for k in passenger[i]))

# if no arc then no transfer or flow
#@constraint(pTran, mustSupply[i in N, j in destination[i], k in passenger[i], t in T; i<j],  trnsfer[i,j,k,t] + trnsfer[j,i,k,t] <= (d[i,j,t]+d[j,i,t])*rhs[k,k])

#@constraint(pTran, mustSupply2[i in N, j in destination[i], k in passenger[i], t in T; i<j], x[i,j,k,t] + x[j,i,k,t] <= (d[i,j,t]+d[j,i,t])*rhs[k,k])

#@constraint(pTran, mustSupply3[i in N, j in destination[i], k in passenger[i], t in T], trnsfer[i,j,k,t] <= x[i,j,k,t])

# passenger flow cost between two stops - objective expression
@expression(pTran, costobj1[i in N, j in destination[i], k in passenger[i], t in T],
                            pth[i,j]*x[i,j,k,t])

# transfer cost
@expression(pTran, costobj2[i in N, j in destination[i], k in passenger[i], t in T],
                            params.penalty*trnsfer[i,j,k,t])

#= line flow cost between two stops - objective expression
@expression(pTran, costobj3[i in N, j in destination[i], t in T],
                            sum(pth[i,j]*d[i,j,t]))
=#

#@expression(pTran, costobj4[t in T],
#                            1000*f[t])

# take preset routes from a file
#= imports a route set from a solution file then fix them
route = CSV.read("/Users/akrm/Dropbox/Julia/transitComm/encapsulated/output/routes.csv", header = true, DataFrame)
for i in 1:size(route,1)
    set_start_value(d[route[i,1], route[i,2], route[i,3]] ,1)
end
=#

# objective function
@objective(pTran, Min, sum(costobj1)
                    #    +  sum(costobj2)
                        )

optimize!(pTran)

#= imports a route set from a solution file then fix them
route = CSV.read("/Users/akrm/Dropbox/Julia/transitComm/encapsulated/output/routes.csv", header = true)
for i in 1:size(route,1)
    set_start_value(d[route[i,1], route[i,2], route[i,3]], 1)
end
=#

#= imports a route set file and warmstart with this solution
route = CSV.read("/Users/akrm/Dropbox/Julia/transitdata/mandl/routesForMandl4.csv", header = false)

for i in 1:size(route,1), j in 1:size(route,2)-1
    if !ismissing(route[i,j+1])
        set_start_value(d[route[i,j], route[i,j+1], i], 1)
    end
end
=#

# Fixing the starting and the ending terminals of the routes
#=
fix(d[16,1,1], 1)
fix(d[11,17,1], 1, force=true)
fix(d[16,9,2], 1, force=true)
fix(d[12,17,2], 1, force=true)
fix(d[16,10,3], 1, force=true)
fix(d[12,17,3], 1, force=true)
fix(d[16,12,4], 1, force=true)
fix(d[1,17,4], 1, force=true)
=#

#=
fix(d[16,1,1], 1)
fix(d[10,17,1], 1, force=true)
fix(d[16,2,2], 1, force=true)
fix(d[11,17,2], 1, force=true)
fix(d[16,5,3], 1, force=true)
fix(d[9,17,3], 1, force=true)
fix(d[16,7,4], 1, force=true)
fix(d[8,17,4], 1, force=true)
=#


#= warm start with our solution
set_start_value(d[1,2,1], 1)
set_start_value(d[2,5,1], 1)
set_start_value(d[5,4,1], 1)
set_start_value(d[4,6,1], 1)
set_start_value(d[6,8,1], 1)
set_start_value(d[8,10,1], 1)
set_start_value(d[10,11,1], 1)

set_start_value(d[9,15,2], 1)
set_start_value(d[15,8,2], 1)
set_start_value(d[8,10,2], 1)
set_start_value(d[10,14,2], 1)
set_start_value(d[14,13,2], 1)
set_start_value(d[13,11,2], 1)
set_start_value(d[11,12,2], 1)

set_start_value(d[10,7,3], 1)
set_start_value(d[7,15,3], 1)
set_start_value(d[15,6,3], 1)
set_start_value(d[6,3,3], 1)
set_start_value(d[3,2,3], 1)
set_start_value(d[2,4,3], 1)
set_start_value(d[4,12,3], 1)

set_start_value(d[12,11,4], 1)
set_start_value(d[11,10,4], 1)
set_start_value(d[10,8,4], 1)
set_start_value(d[8,6,4], 1)
set_start_value(d[6,3,4], 1)
set_start_value(d[3,2,4], 1)
set_start_value(d[2,1,4], 1)
=#




#= objective with only distance solution, transfer is limited, obj 117 - 170000
set_start_value(d[2,4,1],1)
set_start_value(d[3,2,1],1)
set_start_value(d[4,5,1],1)
set_start_value(d[6,3,1],1)
set_start_value(d[7,15,1],1)
set_start_value(d[8,6,1],1)
set_start_value(d[15,8,1],1)
set_start_value(d[6,15,2],1)
set_start_value(d[7,10,2],1)
set_start_value(d[8,6,2],1)
set_start_value(d[10,11,2],1)
set_start_value(d[11,13,2],1)
set_start_value(d[13,14,2],1)
set_start_value(d[15,7,2],1)
set_start_value(d[4,5,3],1)
set_start_value(d[6,4,3],1)
set_start_value(d[7,10,3],1)
set_start_value(d[8,6,3],1)
set_start_value(d[9,15,3],1)
set_start_value(d[10,8,3],1)
set_start_value(d[15,7,3],1)
set_start_value(d[1,2,4],1)
set_start_value(d[2,3,4],1)
set_start_value(d[3,6,4],1)
set_start_value(d[6,8,4],1)
set_start_value(d[8,10,4],1)
set_start_value(d[10,11,4],1)
set_start_value(d[11,12,4],1)
=#

#= objective with only distance solution, transfer is limited, obj 156 - 173250
set_start_value(d[2,5,1],1)
set_start_value(d[4,6,1],1)
set_start_value(d[5,4,1],1)
set_start_value(d[6,8,1],1)
set_start_value(d[8,10,1],1)
set_start_value(d[10,11,1],1)
set_start_value(d[11,13,1],1)

set_start_value(d[4,6,2],1)
set_start_value(d[6,15,2],1)
set_start_value(d[8,10,2],1)
set_start_value(d[10,11,2],1)
set_start_value(d[11,12,2],1)
set_start_value(d[12,4,2],1)
set_start_value(d[15,9,2],1)

set_start_value(d[7,15,3],1)
set_start_value(d[10,7,3],1)
set_start_value(d[11,13,3],1)
set_start_value(d[12,11,3],1)
set_start_value(d[13,14,3],1)
set_start_value(d[14,10,3],1)
set_start_value(d[15,9,3],1)

set_start_value(d[1,2,4],1)
set_start_value(d[2,3,4],1)
set_start_value(d[3,6,4],1)
set_start_value(d[6,8,4],1)
set_start_value(d[7,15,4],1)
set_start_value(d[8,10,4],1)
set_start_value(d[10,7,4],1)
=#


#= warm start with a solution
set_start_value(d[1,2,1], 1)
set_start_value(d[2,5,1], 1)
set_start_value(d[5,4,1], 1)
set_start_value(d[4,6,1], 1)
set_start_value(d[6,8,1], 1)
set_start_value(d[8,10,1], 1)
set_start_value(d[10,7,1], 1)

set_start_value(d[1,2,2], 1)
set_start_value(d[2,4,2], 1)
set_start_value(d[4,6,2], 1)
set_start_value(d[6,15,2], 1)
set_start_value(d[15,7,2], 1)
set_start_value(d[7,10,2], 1)
set_start_value(d[10,14,2], 1)

set_start_value(d[6,4,3], 1)
set_start_value(d[4,12,3], 1)
set_start_value(d[12,11,3], 1)
set_start_value(d[11,10,3], 1)
set_start_value(d[10,8,3], 1)
set_start_value(d[8,15,3], 1)
set_start_value(d[15,9,3], 1)

set_start_value(d[1,2,4], 1)
set_start_value(d[2,3,4], 1)
set_start_value(d[3,6,4], 1)
set_start_value(d[6,8,4], 1)
set_start_value(d[8,10,4], 1)
set_start_value(d[10,13,4], 1)
set_start_value(d[13,11,4], 1)
=#

#= a good solution
set_start_value(d[7,10,1], 1)
set_start_value(d[10,13,1], 1)
set_start_value(d[13,11,1], 1)
set_start_value(d[11,12,1], 1)
set_start_value(d[12,4,1], 1)
set_start_value(d[4,2,1], 1)
set_start_value(d[2,1,1], 1)

set_start_value(d[7,15,2], 1)
set_start_value(d[15,8,2], 1)
set_start_value(d[8,6,2], 1)
set_start_value(d[6,4,2], 1)
set_start_value(d[4,5,2], 1)
set_start_value(d[5,2,2], 1)
set_start_value(d[2,1,2], 1)

set_start_value(d[1,2,3], 1)
set_start_value(d[2,3,3], 1)
set_start_value(d[3,6,3], 1)
set_start_value(d[6,15,3], 1)
set_start_value(d[15,7,3], 1)
set_start_value(d[7,10,3], 1)
set_start_value(d[10,14,3], 1)

set_start_value(d[9,15,4], 1)
set_start_value(d[15,8,4], 1)
set_start_value(d[8,6,4], 1)
set_start_value(d[6,3,4], 1)
set_start_value(d[3,2,4], 1)
set_start_value(d[2,4,4], 1)
set_start_value(d[4,12,4], 1)

set_start_value(d[5,4,5], 1)
set_start_value(d[4,12,5], 1)
set_start_value(d[12,11,5], 1)
set_start_value(d[11,10,5], 1)
set_start_value(d[10,7,5], 1)
set_start_value(d[7,15,5], 1)
set_start_value(d[15,9,5], 1)

set_start_value(d[9,15,6], 1)
set_start_value(d[15,6,6], 1)
set_start_value(d[6,4,6], 1)

set_start_value(d[1,2,7], 1)
set_start_value(d[2,3,7], 1)
set_start_value(d[3,6,7], 1)
set_start_value(d[6,8,7], 1)
set_start_value(d[8,10,7], 1)
set_start_value(d[10,11,7], 1)
set_start_value(d[11,13,7], 1)

set_start_value(d[14,13,8], 1)
set_start_value(d[13,11,8], 1)
set_start_value(d[11,10,8], 1)
set_start_value(d[10,8,8], 1)
set_start_value(d[8,6,8], 1)
set_start_value(d[6,4,8], 1)
set_start_value(d[4,5,8], 1)
=#


#= warm start with basic Benders solution
set_start_value(d[1,2,1], 1)
set_start_value(d[2,5,1], 1)
set_start_value(d[5,4,1], 1)
set_start_value(d[4,12,1], 1)
set_start_value(d[12,11,1], 1)
set_start_value(d[11,10,1], 1)
set_start_value(d[10,8,1], 1)

set_start_value(d[12,4,2], 1)
set_start_value(d[4,5,2], 1)
set_start_value(d[5,2,2], 1)
set_start_value(d[2,3,2], 1)
set_start_value(d[3,6,2], 1)
set_start_value(d[6,15,2], 1)
set_start_value(d[15,9,2], 1)

set_start_value(d[11,13,3], 1)
set_start_value(d[13,14,3], 1)
set_start_value(d[14,10,3], 1)
set_start_value(d[10,8,3], 1)
set_start_value(d[8,15,3], 1)
set_start_value(d[15,6,3], 1)
set_start_value(d[6,4,3], 1)

set_start_value(d[15,8,4], 1)
set_start_value(d[8,10,4], 1)
set_start_value(d[10,14,4], 1)
set_start_value(d[14,13,4], 1)
set_start_value(d[13,11,4], 1)
set_start_value(d[11,12,4], 1)
set_start_value(d[12,4,4], 1)
=#
