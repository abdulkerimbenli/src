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
@variable(subProb, mu[i in N, j in destination[i], t in T], lower_bound=0)

# dual variable for primal problem pTran's transfer constraints
@variable(subProb, teta[i in N, j in destination[i], k in passenger[i], t in T], lower_bound=0)

# dual variable for primal problem pTran's maximum itinerary constraint - free variable
@variable(subProb, beta[t in T])

@variable(subProb, dual_artificial[i in N, j in destination[i], t in T])

@variable(subProb, d[i in N, j in destination[i], t in T], Bin)

# constraint for primal problem pTran's passenger flow variable
@constraint(subProb, dualFlow[i in N, j in destination[i], k in passenger[i], t in T; i==k],
                            lambda[i,k] - lambda[j,k]
                            + sum(teta[j,r,k,t] for r in destination[j] if r!=i)
                            - mu[i,j,t] <= pth[i,j])

# constraint for primal problem pTran's passenger flow variable
@constraint(subProb, dualFlow2[i in N, j in destination[i], k in passenger[i], t in T; i!=k],
                            lambda[i,k] - lambda[j,k]
                            - teta[i,j,k,t]
                            + sum(teta[j,r,k,t] for r in destination[j] if r!=i)
                            - mu[i,j,t] <= pth[i,j])

# constraint for primal problem pTran's transfer variable
@constraint(subProb, dualTransfer[i in N, j in destination[i], k in passenger[i], t in T],
                            teta[i,j,k,t] <= params.penalty)

# constraint for primal problem pTran's line flow variable (d[i,j,t])
@constraint(subProb, dualLink[i in N, j in destination[i], t in T],
                            params.bigNumber*(mu[i,j,t] + mu[j,i,t])
                            +dual_artificial[i,j,t] == 0)

# objective term for  primal problem pTran's passenger flow conservation constraint
ex1 = @expression(subProb, dualObj1, sum(rhs[i,k]*lambda[i,k] for i in N, k in passenger[i]))

# objective term for primal problem pTran's route itinerary constraint
#ex2 = @expression(subProb, dualObj2, sum((params.itinerary-1)*beta[t] for t in T))



# test the optimal solution
# take preset routes from a file
# imports a route set from a solution file then fix them
d_current = zeros((params.nNodes+2, params.nNodes+2, params.nLine))

route = CSV.read("/Users/akrm/Dropbox/Julia/transitComm/encapsulated/output/routes.csv", header = true, DataFrame)
for i in 1:size(route,1)
    d_current[route[i,1], route[i,2], route[i,3]] = 1
end
#


# objective term for primal Problme pTran's capacity-link constraints

ex3 = @expression(subProb, dualObj3, sum(dual_artificial[i,j,t]*d_current[i,j,t]
                                          for i in N, j in destination[i], t in T))

@objective(subProb, Max, ex1 + ex2 + ex3)

optimize!(subProb)
