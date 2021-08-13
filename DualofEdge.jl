## DUAL PROBLEM
#

# import data
flowq = ODistance(true, false)
pth = ODistance(false, true)
sup, rhs = rightHandSide()
params = problemParameters()
N, NT, T, destination, passenger, lineDict = requiredDictionaries()

# we use direct_model function since we will deal with the backend of the model
# in the backend we will retrieve the extreme rays of unbounded subproblem cases
subProb = direct_model(Gurobi.Optimizer())

# we need info about the unbounded subproblem
set_optimizer_attribute(subProb,"InfUnbdInfo",1)
set_optimizer_attribute(subProb,"Presolve",0)
set_optimizer_attribute(subProb,"Heuristics",0.0)
set_optimizer_attribute(subProb,"Cuts",0)

# dual variable for primal problem pTran's passenger flow conservation constraints - free variable
@variable(subProb, lambda[i in N, k in passenger[i]])

# dual variable for primal problem pTran's linking constraints with d_ijt
@variable(subProb, mu[i in N, j in destination[i], t in T; i<j], upper_bound=0)

# dual variable for primal problem pTran's transfer constraints
@variable(subProb, teta[i in N, j in destination[i], k in passenger[i], t in T; i!=k], upper_bound=0)

# dual variable for primal problem pTran's maximum itinerary constraint - free variable
@variable(subProb, beta[t in T])

# constraint for primal problem pTran's passenger flow variable
@constraint(subProb, dualFlow[i in N, j in destination[i], k in passenger[i], t in T; j!=k && i==k && i<j],
                            lambda[i,k] - lambda[j,k]
                            - sum(teta[j,r,k,t] for r in destination[j] if r!=i)
                            + mu[i,j,t] <= pth[i,j])

@constraint(subProb, dualFlowA[i in N, j in destination[i], k in passenger[i], t in T; j!=k && i==k && j<i],
                            lambda[i,k] - lambda[j,k]
                            - sum(teta[j,r,k,t] for r in destination[j] if r!=i)
                            + mu[j,i,t] <= pth[i,j])

# constraint for primal problem pTran's passenger flow variable
@constraint(subProb, dualFlow2[i in N, j in destination[i], k in passenger[i], t in T; j!=k && i!=k && i<j],
                            lambda[i,k] - lambda[j,k]
                            + teta[i,j,k,t]
                            - sum(teta[j,r,k,t] for r in destination[j] if r!=i)
                            + mu[i,j,t] <= pth[i,j])

@constraint(subProb, dualFlow2A[i in N, j in destination[i], k in passenger[i], t in T; j!=k && i!=k && j<i],
                            lambda[i,k] - lambda[j,k]
                            + teta[i,j,k,t]
                            - sum(teta[j,r,k,t] for r in destination[j] if r!=i)
                            + mu[j,i,t] <= pth[i,j])

@constraint(subProb, dualFlow3[i in N, j in destination[i], k in passenger[i], t in T; j==k && i!=k && i<j],
                            lambda[i,k]
                            + teta[i,j,k,t]
                            + mu[i,j,t] <= pth[i,j])

@constraint(subProb, dualFlow3A[i in N, j in destination[i], k in passenger[i], t in T; j==k && i!=k && j<i],
                            lambda[i,k]
                            + teta[i,j,k,t]
                            + mu[j,i,t] <= pth[i,j])

# constraint for primal problem pTran's line flow variable (d[i,j,t])
@constraint(subProb, dualLink[i in N, j in destination[i], t in T; i<j],
                            - params.bigNumber*(mu[i,j,t])
                            + beta[t] <= 0)

# constraint for primal problem pTran's transfer variable
@constraint(subProb, dualTransfer[i in N, j in destination[i], k in passenger[i], t in T; i!=k],
                            - teta[i,j,k,t] <= params.penalty)


# objective term for  primal problem pTran's passenger flow conservation constraint
ex1 = @expression(subProb, dualObj1, sum(rhs[i,k]*lambda[i,k] for i in N, k in passenger[i]))

# objective term for primal problem pTran's route itinerary constraint
ex2 = @expression(subProb, dualObj2, sum((params.itinerary-1)*beta[t] for t in T))


d_current = zeros((params.nNodes, params.nNodes, params.nLine))

# test the optimal solution
d_current[1,2,1] = 1
d_current[2,5,1] = 1
d_current[4,5,1] = 1
d_current[4,6,1] = 1
d_current[6,8,1] = 1
d_current[8,10,1] = 1
#d_current[10,11,1] = 1

d_current[9,15,2] = 1
d_current[8,15,2] = 1
d_current[8,10,2] = 1
d_current[10,14,2] = 1
d_current[13,14,2] = 1
d_current[11,13,2] = 1
d_current[11,12,2] = 1

d_current[7,10,3] = 1
d_current[7,15,3] = 1
d_current[6,15,3] = 1
d_current[3,6,3] = 1
d_current[2,3,3] = 1
d_current[2,4,3] = 1
d_current[4,12,3] = 1

d_current[11,12,4] = 1
d_current[10,11,4] = 1
d_current[8,10,4] = 1
d_current[6,8,4] = 1
d_current[3,6,4] = 1
d_current[2,3,4] = 1
d_current[1,2,4] = 1
#

#=
route = CSV.read("/Users/akrm/Dropbox/Julia/transitComm/encapsulated/output/routes.csv", header = true, DataFrame)
for i in 1:size(route,1)
    d_current[route[i,1], route[i,2], route[i,3]] = 1
end
=#

#objective term for primal Problme pTran's capacity-link constraints

ex3 = @expression(subProb, dualObj3, sum((params.bigNumber*mu[i,j,t] - beta[t])*d_current[i,j,t]
                                          for i in N, j in destination[i], t in T if i<j))

@objective(subProb, Max, ex1 + ex2 + ex3)


optimize!(subProb)


## PRIMAL PROBLEM
#
#

# if there are number of exact number of itineraries then the model gives arbitrary numbers for d values and
# the optimal value cannot be found
# so the maximum number of itinerary must be fixed
# AND BE CAREFUL HERE, the feasibility of arcs ( d variable) is not dealt here, so results can be infeasible in the original problem

primal = Model(optimizer_with_attributes(Gurobi.Optimizer, "Presolve" => 0))


# passenger flow decision variable
@variable(primal, x[i in N, j in destination[i], k in passenger[i], t in T], lower_bound=0)

@variable(primal, d[i in N, j in destination[i], t in T; i<j])

# transfer flow decision variable
@variable(primal, trnsfer[i in N, j in destination[i], k in passenger[i], t in T], lower_bound=0)


# PASSENGER FLOW CONSERVATION CONSTRAINTS:
# supply of passenger flows
@constraint(primal, flowb[i in N],
                            sum(x[i,j,i,t] for j in destination[i], t in T) == rhs[i,i])

# passenger flows leaves lines at the desired stop, otherwise just pass over that stops
@constraint(primal, flowb2[i in N, k in passenger[i]; k!=i],
                            + sum(x[i,j,k,t] for j in destination[i], t in T)
                            - sum(x[j,i,k,t] for j in destination[i], t in T) == rhs[i,k])

# passenger arrives a stop from the previous stop of the line or
# transfer from another line (except the stop's own passengers)
@constraint(primal, transfer[i in N, j in destination[i], k in passenger[i], t in T; (k!=i)],
                            x[i,j,k,t] <= sum(x[r,i,k,t] for r in destination[i] if r!=j)
                                                + trnsfer[i,j,k,t])

# a passenger can through one stop to another iff there is line flow over there
# in other words the segment of a line is defined
@constraint(primal, capacity[i in N, j in destination[i], k in passenger[i], t in T; i<j],
                            x[i,j,k,t] <= params.bigNumber*d[i,j,t])

@constraint(primal, capacity2[i in N, j in destination[i], k in passenger[i], t in T; j<i],
                            x[i,j,k,t] <= params.bigNumber*d[j,i,t])

#

# passenger flow cost between two stops - objective expression
@expression(primal, costobj1[i in N, j in destination[i], k in passenger[i], t in T],
                            pth[i,j]*x[i,j,k,t])

# transfer cost
@expression(primal, costobj2[i in N, j in destination[i], k in passenger[i], t in T],
                            params.penalty*trnsfer[i,j,k,t])

@constraint(primal, link4[t in T],
                        sum(d[i,j,t] for i in N,  j in destination[i] if i<j)
                        <= params.itinerary-1)

# warm start with our solution
fix(d[1,2,1], 1; force = true)
fix(d[2,5,1], 1; force = true)
fix(d[4,5,1], 1; force = true)
fix(d[4,6,1], 1; force = true)
fix(d[6,8,1], 1; force = true)
fix(d[8,10,1], 1; force = true)
fix(d[10,11,1], 1; force = true)

fix(d[9,15,2], 1; force = true)
fix(d[8,15,2], 1; force = true)
fix(d[8,10,2], 1; force = true)
fix(d[10,14,2], 1; force = true)
fix(d[13,14,2], 1; force = true)
fix(d[11,13,2], 1; force = true)
fix(d[11,12,2], 1; force = true)

fix(d[7,10,3], 1; force = true)
fix(d[7,15,3], 1; force = true)
fix(d[6,15,3], 1; force = true)
fix(d[3,6,3], 1; force = true)
fix(d[2,3,3], 1; force = true)
fix(d[2,4,3], 1; force = true)
fix(d[4,12,3], 1; force = true)

fix(d[11,12,4], 1; force = true)
fix(d[10,11,4], 1; force = true)
fix(d[8,10,4], 1; force = true)
fix(d[6,8,4], 1; force = true)
fix(d[3,6,4], 1; force = true)
fix(d[2,3,4], 1; force = true)
fix(d[1,2,4], 1; force = true)
#


# objective function
@objective(primal, Min, sum(costobj1) + sum(costobj2))
#
optimize!(primal)
#


# TAKE THE DULA OF THE PRIMAL MODEL
using Dualization

dual_model = dualize(primal; dual_names = DualNames("dual", ""))

set_optimizer(dual_model, Gurobi.Optimizer)

optimize!(dual_model)

open("dualEdge.txt", "w") do f
   println(f, dual_model)
end

#
