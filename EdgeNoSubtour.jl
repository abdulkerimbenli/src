"""
an experimental model with the edge formulation
"""

flowq = ODistance(true,false)
pth = ODistance(false,true)
sup, rhs = rightHandSide()
params = problemParameters()
N, NT, T, destination, passenger,lineDict = requiredDictionaries()
# mathematical model

pTran = Model(optimizer_with_attributes(Gurobi.Optimizer,
                #"Presolve" => 0,
                #"MIPGap" => 0.0001,
                #"MIPFocus" => 1,
                #"Heuristics" => 0.95,
                #"Cuts" => 0,
                #"ImproveStartTime" => 1000,
                #"Threads" => 8,
                #"TimeLimit" => 100
                ))

# add norelheur heuristics in which heusristics are used before relaxation for the cases relaxation takes time
#set_optimizer_attribute(pTran,"NoRelHeurTime",200)
#

# passenger flow decision variable
@variable(pTran, x[i in N, j in destination[i], k in passenger[i], t in T], lower_bound=0)

# transfer flow decision variable
@variable(pTran, trnsfer[i in N, j in destination[i], k in passenger[i], t in T], lower_bound=0)

# line flow decision variable
@variable(pTran, d[i in NT, j in lineDict[i], t in T; i<j || i>length(N)], Bin)

# adding variable frequency of lines
@variable(pTran, f[t in T], lower_bound=0.0)

# starting dummny node can only connect one node in each route
@constraint(pTran, link1[i in length(N)+1, t in T],
                            sum(d[i,j,t] for j in lineDict[i]) == 1)

# sinking dummny node can only connect one node in each route
@constraint(pTran, link2[i in length(N)+2, t in T],
                            sum(d[j,i,t] for j in N) == 1)

# create a binary variable to control the right-hand side of the edge constraint
@variable(pTran, z[i in N, t in T], Bin)

# since we use edge formulation a node (except dummies) can connect two nodes
@constraint(pTran, link3[i in N, t in T],
                        sum(d[i,j,t] for j in lineDict[i] if i<j)
                        + sum(d[j,i,t] for j in lineDict[i] if j<i)
                        + sum(d[length(N)+1, i, t]) == 2*z[i,t])

# in a route there are only a number of arcs including links which connects dummies to the network
# this way we try to tackle the subtour problem
@constraint(pTran, link4[t in T],
                        sum(d[i,j,t] for i in NT,  j in lineDict[i] if (i<j || i in length(N)+1))
                        <= params.itinerary+1)

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

#= put range for number of passengers
@constraint(pTran, numTransfer, sum(trnsfer[i,j,k,t]
                                    for i in N, j in destination[i], k in passenger[i], t in T)
                                        <= 0.081*sum(sup))
=#

# a passenger can through one stop to another iff there is line flow over there
# in other words the segment of a line is defined
@constraint(pTran, capacity[i in N, j in destination[i], t in T; i<j],
                            sum(x[i,j,k,t] for k in passenger[i])
                                <=params.bigNumber*params.cap[t]*d[i,j,t])

@constraint(pTran, capacity2[i in N, j in destination[i], t in T; j<i],
                            sum(x[i,j,k,t] for k in passenger[i])
                                <=params.bigNumber*params.cap[t]*d[j,i,t])

# frequency determination constraint
@constraint(pTran, freq[i in N, j in destination[i], t in T],
                           sum(x[i,j,k,t] for k in passenger[i]) + sum(x[j,i,k,t] for k in passenger[j])
                                <= f[t]*params.cap[t])

# passenger flow cost between two stops - objective expression
@expression(pTran, costobj1[i in N, j in destination[i], k in passenger[i], t in T],
                            pth[i,j]*x[i,j,k,t])

# transfer cost
@expression(pTran, costobj2[i in N, j in destination[i], k in passenger[i], t in T],
                            params.penalty*trnsfer[i,j,k,t])

# objective function
@objective(pTran, Min,  sum(costobj1)
                        +  sum(costobj2)
                        #+ sum(costobj3)
                        )

#= imports a route set from a solution file then fix them
route = CSV.read("/Users/akrm/Dropbox/Julia/transitComm/encapsulated/output/routes.csv", header = true, DataFrame)
for i in 1:size(route,1)
    set_start_value(d[route[i,1], route[i,2], route[i,3]], 1)
    #println(route[i,1], ",", route[i,2], ",", route[i,3])
end
=#

#=
d_current = zeros((params.nNodes, params.nNodes, params.nLine))

#route = CSV.read("/Users/akrm/Dropbox/Julia/transitComm/encapsulated/output/routes.csv", header = true, DataFrame)
for i in 1:size(route,1)
    d_current[route[i,1], route[i,2], route[i,3]] = 1
end
#

for i in N, j in destination[i], t in T
    if i<j
        fix(d[i,j,t], d_current[i,j,t])
    end
end

=#

optimize!(pTran)

#= warm start with our solution
set_start_value(d[1,2,1], 1)
set_start_value(d[2,5,1], 1)
set_start_value(d[4,5,1], 1)
set_start_value(d[4,6,1], 1)
set_start_value(d[6,8,1], 1)
set_start_value(d[8,10,1], 1)
set_start_value(d[10,11,1], 1)

set_start_value(d[9,15,2], 1)
set_start_value(d[8,15,2], 1)
set_start_value(d[8,10,2], 1)
set_start_value(d[10,14,2], 1)
set_start_value(d[13,14,2], 1)
set_start_value(d[11,13,2], 1)
set_start_value(d[11,12,2], 1)

set_start_value(d[7,10,3], 1)
set_start_value(d[7,15,3], 1)
set_start_value(d[6,15,3], 1)
set_start_value(d[3,6,3], 1)
set_start_value(d[2,3,3], 1)
set_start_value(d[2,4,3], 1)
set_start_value(d[4,12,3], 1)

set_start_value(d[11,12,4], 1)
set_start_value(d[10,11,4], 1)
set_start_value(d[8,10,4], 1)
set_start_value(d[6,8,4], 1)
set_start_value(d[3,6,4], 1)
set_start_value(d[2,3,4], 1)
set_start_value(d[1,2,4], 1)
=#
