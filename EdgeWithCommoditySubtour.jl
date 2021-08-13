"""
we try to eliminate the subtours in the edge formulation
by introducing the line commodities which will traverse only self own lines

"""

flowq = ODistance(true,false)
pth = ODistance(false,true)
sup, rhs = rightHandSide()
params = problemParameters()
N, NT, T, destination, passenger,lineDict = requiredDictionaries()

# mathematical model
pTran = Model(optimizer_with_attributes(Gurobi.Optimizer,
                "Presolve" => 0,
                #"MIPGap" => 0.0001,
                #"MIPFocus" => 1,
                #"Heuristics" => 0.95,
                #"Cuts" => 0,
                #"ImproveStartTime" => 1000,
                #"Threads" => 8,
                #"TimeLimit" => 100
                ))

# add norelheur heuristics in which heusristics are used before relaxation for the cases relaxation takes time
#set_optimizer_attribute(pTran,"NoRelHeurTime",10)
#

# the new commodity set including line commodities
K = [N; collect(params.nNodes+3:params.nNodes+3+params.nLine-1)]

# passenger flow decision variable
# a line commodity can only traverse its route
@variable(pTran, x[i in NT, j in lineDict[i], k in K, t in T; k<=params.nNodes || k==t+params.nNodes+2], lower_bound=0)

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

# flow of line commodities from the super source node to super sink node
@constraint(pTran, flowLine[i in params.nNodes+1, k in params.nNodes+3:params.nNodes+3+params.nLine-1, t in T; k==t+params.nNodes+2],
                            sum(x[i,j,k,t] for j in lineDict[i]) == 1)

@constraint(pTran, flowLine2[i in N, k in params.nNodes+3:params.nNodes+3+params.nLine-1, t in T; k==t+params.nNodes+2],
                            +sum(x[i,j,k,t] for j in lineDict[i])
                            -sum(x[j,i,k,t] for j in destination[i])
                            -x[params.nNodes+1, i, k,t] == 0)

@constraint(pTran, flowLine3[k in params.nNodes+3:params.nNodes+3+params.nLine-1, t in T; k==t+params.nNodes+2],
                            -sum(x[i,params.nNodes+2,k,t] for i in N) == -1)

#= THEY CAN BE UNLOCKED IT IT IS NECESSARY
@constraint(pTran, linePath[i in NT, j in lineDict[i], k in params.nNodes+3:params.nNodes+3+params.nLine-1, t in T;
                                    (i<j || i>length(N)) && k==t+params.nNodes+2],
                            x[i,j,k,t] <= d[i,j,t])

@constraint(pTran, linePath2[i in N, j in lineDict[i], k in params.nNodes+3:params.nNodes+3+params.nLine-1, t in T;
                                    j<i && k==t+params.nNodes+2],
                            x[i,j,k,t] <= d[j,i,t])
=#

# ensure for subtour elimination,
#= route commodity only traverse if the path is direct i.e. without subtours
@constraint(pTran, finalTouch[i in N, j in destination[i], k in params.nNodes+3:params.nNodes+3+params.nLine-1, t in T;
                                    i<j && k==t+params.nNodes+2],
                            x[i,j,k,t] + x[j,i,k,t] == d[i,j,t])

@constraint(pTran, finalTouch[i in N, j in destination[i], k in params.nNodes+3:params.nNodes+3+params.nLine-1, t in T;
                                    k==t+params.nNodes+2],
                            x[i,j,k,t] == (z[i,t] + z[j,t])*0.5)

@constraint(pTran, finalTouch[i in N, j in destination[i], k in params.nNodes+3:params.nNodes+3+params.nLine-1, t in T;
                                    i<j && k==t+params.nNodes+2],
                                    d[i,j,t] == x[i,j,k,t] + x[j,i,k,t])

@constraint(pTran, finalTouch2[i in N, j in destination[i], k in params.nNodes+3:params.nNodes+3+params.nLine-1, t in T;
                                    k==t+params.nNodes+2],
                                    x[i,j,k,t] <= 0.5)

=#



# passenger arrives a stop from the previous stop of the line or
# transfer from another line (except the stop's own passengers)
@constraint(pTran, transfer[i in N, j in destination[i], k in passenger[i], t in T; (k!=i)],
                            x[i,j,k,t] <= sum(x[r,i,k,t] for r in destination[i] if r!=j)
                                                + trnsfer[i,j,k,t])

#= put range for number of passengers
@constraint(pTran, numTransfer, sum(trnsfer[i,j,k,t]
                                    for i in N, j in destination[i], k in passenger[i], t in T)
                                        == 0.10*params.bigNumber*params.cap[1])
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

#= warm start with our solution
fix(d[1,2,1], 1)
fix(d[2,5,1], 1)
fix(d[4,5,1], 1)
fix(d[4,6,1], 1)
fix(d[6,8,1], 1)
fix(d[8,10,1], 1)
fix(d[10,11,1], 1)

fix(d[9,15,2], 1)
fix(d[8,15,2], 1)
fix(d[8,10,2], 1)
fix(d[10,14,2], 1)
fix(d[13,14,2], 1)
fix(d[11,13,2], 1)
fix(d[11,12,2], 1)

fix(d[7,10,3], 1)
fix(d[7,15,3], 1)
fix(d[6,15,3], 1)
fix(d[3,6,3], 1)
fix(d[2,3,3], 1)
fix(d[2,4,3], 1)
fix(d[4,12,3], 1)

fix(d[11,12,4], 1)
fix(d[10,11,4], 1)
fix(d[8,10,4], 1)
fix(d[6,8,4], 1)
fix(d[3,6,4], 1)
fix(d[2,3,4], 1)
fix(d[1,2,4], 1)
=#

# objective function
@objective(pTran, Min,  sum(costobj1) + sum(costobj2) )
optimize!(pTran)

#=
using Dualization

dual_model = dualize(pTran; dual_names = DualNames("dual", ""))

set_optimizer(dual_model, Gurobi.Optimizer)

optimize!(dual_model)

print(objective_function(dual_model))
=#
