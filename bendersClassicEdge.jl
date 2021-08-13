"""
we try to eliminate the subtours in the edge formulation
by introducing the line commodities which will traverse self own lines

"""

flowq = ODistance(true, false)
pth = ODistance(false, true)
sup, rhs = rightHandSide()
params = problemParameters()
N, NT, T, destination, passenger, lineDict = requiredDictionaries()
# the new commodity set including line commodities
# K = [N; collect(params.nNodes+3:params.nNodes+3+params.nLine-1)]

master = direct_model(Gurobi.Optimizer())

set_optimizer_attribute(master,"Presolve",0)
set_optimizer_attribute(master,"Heuristics",0.0)
set_optimizer_attribute(master,"Cuts",0)
#set_optimizer_attribute(master, "Threads",1) # sometimes parallel cores find solutions in other nodes and lead to wrong answers
#set_optimizer_attribute(master,"Method",0)

# create a surrogate variable to hold subProblem's objective in the master
@variable(master, p, lower_bound = -1e9)

# line flow decision variable
@variable(master, d[i in NT, j in lineDict[i], t in T; i<j || i>length(N)], Bin)

# starting dummny node can only connect one node in each route
@constraint(master, link1[i in length(N)+1, t in T],
                            sum(d[i,j,t] for j in lineDict[i]) == 1)

# sinking dummny node can only connect one node in each route
@constraint(master, link2[i in length(N)+2, t in T],
                            sum(d[j,i,t] for j in N) == 1)

# create a binary variable to control the right-hand side of the edge constraint
@variable(master, z[i in N, t in T], Bin)

# since we use edge formulation a node (except dummies) can connect two nodes
@constraint(master, link3[i in N, t in T],
                        sum(d[i,j,t] for j in lineDict[i] if i<j)
                        + sum(d[j,i,t] for j in lineDict[i] if j<i)
                        + sum(d[length(N)+1, i, t]) == 2 * z[i,t])

@constraint(master, limitItinerary[t in T],
                        sum(d[i,j,t] for i in N, j in destination[i] if i<j) == params.itinerary-1)

#
@objective(master, Min, p)

#
fix(d[1,2,1], 1)
#fix(d[2,5,1], 1)
#fix(d[4,5,1], 1)
#fix(d[4,6,1], 1)
#fix(d[6,8,1], 1)
#fix(d[8,10,1], 1)
fix(d[7,10,1], 1)

#
fix(d[9,15,2], 1)
#fix(d[8,15,2], 1)
#fix(d[8,10,2], 1)
#fix(d[10,14,2], 1)
#fix(d[13,14,2], 1)
#fix(d[11,13,2], 1)
#fix(d[11,12,2], 1)
#

fix(d[7,10,3], 1)
#fix(d[7,15,3], 1)
#fix(d[6,15,3], 1)
#fix(d[3,6,3], 1)
#fix(d[2,3,3], 1)
#fix(d[2,4,3], 1)
fix(d[4,12,3], 1)

fix(d[11,12,4], 1)
#set_start_value(d[10,11,4], 1)
#set_start_value(d[8,10,4], 1)
#set_start_value(d[6,8,4], 1)
#set_start_value(d[3,6,4], 1)
#set_start_value(d[2,3,4], 1)
#set_start_value(d[1,2,4], 1)
#

global iter_num = 1

iter_num = 1

while true
    println("\n-----------------------")
    println("Iteration number = ", iter_num)
    println("-----------------------\n")

    #print(master)

    optimize!(master)

    t_status = termination_status(master)
    p_status = primal_status(master)

    if p_status == MOI.INFEASIBLE_POINT
        println("The problem is infeasible :-(")
        break
    end

    if p_status == MOI.FEASIBLE_POINT
    master_current = value(p)
    d_current = value.(d)
    d_current = round.(d_current)
    end

    println(
        "Status of the master problem is ", t_status,"\nwith master_current = ",master_current,
        #"\nx_current = ", x_current
    )

    # we use direct_model function since we will deal with the backend of the model
    # in the backend we will retrieve the extreme rays of unbounded subproblem cases
    subProb = direct_model(Gurobi.Optimizer())

    # we need info about the unbounded subproblem
    set_optimizer_attribute(subProb,"InfUnbdInfo",1)
    set_optimizer_attribute(subProb,"Presolve",0)
    #set_optimizer_attribute(subProb,"Heuristics",0.0)
    #set_optimizer_attribute(subProb,"Cuts",0)
    #set_optimizer_attribute(subProb,"Method",0)
    #set_optimizer_attribute(subProb,"FeasibilityTol", 1e-2)

    # dual variable for primal problem pTran's passenger flow conservation constraints - free variable
    @variable(subProb, lambda[i in N, k in passenger[i]])

    # dual variable for primal problem pTran's linking constraints with d_ijt
    @variable(subProb, mu[i in N, j in destination[i], t in T; i<j], upper_bound=0)

    # dual variable for primal problem pTran's transfer constraints
    @variable(subProb, teta[i in N, j in destination[i], k in passenger[i], t in T; i!=k], upper_bound=0)

    # dual variable for primal problem pTran's maximum itinerary constraint - free variable
    @variable(subProb, beta[t in T])

    # dual variable for the 3rd term of the objective function of the subProblem
    @variable(subProb, artificial[i in N, j in destination[i], t in T; i<j])

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
                            + beta[t] + artificial[i,j,t] == 0)

    # constraint for primal problem pTran's transfer variable
    @constraint(subProb, dualTransfer[i in N, j in destination[i], k in passenger[i], t in T; i!=k],
                            - teta[i,j,k,t] <= params.penalty)


    # objective term for  primal problem pTran's passenger flow conservation constraint
    ex1 = @expression(subProb, dualObj1, sum(rhs[i,k]*lambda[i,k] for i in N, k in passenger[i]))

    # objective term for primal problem pTran's route itinerary constraint
    ex2 = @expression(subProb, dualObj2, sum((params.itinerary-1)*beta[t] for t in T))

    ex3 = @expression(subProb, sum(artificial[i,j,t]*d_current[i,j,t]
                                            for i in N, j in destination[i], t in T if i<j))

    @objective(subProb, Max, ex1 + ex2 + ex3)

    #print(value.(d))
    #print(subProb)

    optimize!(subProb)

    t_status_sub = termination_status(subProb)
    p_status_sub = primal_status(subProb)


    sub_current = objective_value(subProb)

    if p_status_sub == MOI.FEASIBLE_POINT && sub_current == master_current # we are done
        println("\n################################################")
        println("Optimal solution of the original problem found")
        println("The optimal objective value p is ", master_current)
        break
    end

    if p_status_sub == MOI.FEASIBLE_POINT && sub_current > master_current
        println(
            "\nThere is a suboptimal vertex, add the corresponding constraint",
        )

        #construct terms in the optimality cut
        ex4 = @expression(master, sum(rhs[i,k]*value(lambda[i,k]) for i in N, k in passenger[i]))
        ex5 = @expression(master, sum((params.itinerary-1)*value(beta[t]) for t in T))
        ex6 = @expression(master, sum(value(artificial[i,j,t])*d[i,j,t]
                                                  for i in N, j in destination[i], t in T if i<j))

        @constraint(master, p >= ex4 + ex5 + ex6 )
    end

    if t_status_sub == MOI.DUAL_INFEASIBLE
        println(
            "\nThere is an  extreme ray, adding the corresponding constraint",
        )

        # in the backend we will retrieve the extreme rays of unbounded subproblem cases
        grb = backend(subProb)

        # to store extreme rays corresponding the subproblem
        global lambdaF = zeros((params.nNodes,params.nNodes))
        global betaF = zeros(params.nLine)
        global artificialF = zeros((params.nNodes,params.nNodes,params.nLine))

        for i in N, j in destination[i], k in passenger[i], t in T
        if i<j
          r = Ref{Cdouble}()
          a = lambda[i,k].index.value
          Gurobi.GRBgetdblattrelement(grb, GRB_DBL_ATTR_UNBDRAY,a,r)
          lambdaF[i,k] = r.x

          r = Ref{Cdouble}()
          a = beta[t].index.value
          Gurobi.GRBgetdblattrelement(grb, GRB_DBL_ATTR_UNBDRAY,a,r)
          betaF[t] = r.x

          r = Ref{Cdouble}()
          a = artificial[i,j,t].index.value
          Gurobi.GRBgetdblattrelement(grb, GRB_DBL_ATTR_UNBDRAY,a,r)
          artificialF[i,j,t] = r.x
        end
        end

        ex7 = @expression(master, sum(rhs[i,k]*lambdaF[i,k] for i in N, k in passenger[i]))
        ex8 = @expression(master, sum((params.itinerary-1)*betaF[t] for t in T))
        ex9 = @expression(master, sum(artificialF[i,j,t]*d[i,j,t]
                                                for i in N, j in destination[i], t in T if i<j))

        @constraint(master, ex7 + ex8 + ex9 <= 0 )
        println("\n feasibility cut is added to the master")
    end


    global iter_num += 1
end


#=
d_current = zeros((params.nNodes, params.nNodes, params.nLine))

# warm start with our solution
d_current[1,2,1] = 1
d_current[2,5,1] = 1
d_current[5,4,1] = 1
d_current[4,6,1] = 1
d_current[6,8,1] = 1
d_current[8,10,1] = 1
d_current[10,11,1] = 1

d_current[9,15,2] = 1
d_current[15,8,2] = 1
d_current[8,10,2] = 1
d_current[10,14,2] = 1
d_current[14,13,2] = 1
d_current[13,11,2] = 1
d_current[11,12,2] = 1

d_current[10,7,3] = 1
d_current[7,15,3] = 1
d_current[15,6,3] = 1
d_current[6,3,3] = 1
d_current[3,2,3] = 1
d_current[2,4,3] = 1
d_current[4,12,3] = 1

d_current[12,11,4] = 1
d_current[11,10,4] = 1
d_current[10,8,4] = 1
d_current[8,6,4] = 1
d_current[6,3,4] = 1
d_current[3,2,4] = 1
d_current[2,1,4] = 1
=#
