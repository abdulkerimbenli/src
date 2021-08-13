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

#set_optimizer_attribute(master, "Presolve", 2)
#set_optimizer_attribute(master,"Heuristics",0.0)
#set_optimizer_attribute(master,"Cuts",0)
#set_optimizer_attribute(master,"Method",0)

# create a surrogate variable o hold subProblem's objective in the master
@variable(master, p, lower_bound = -1e7)

# line flow decision variable
@variable(master, d[i in NT, j in lineDict[i], t in T], Bin)

# adding MTZ subtour elimination variable
@variable(master, u[i in N, t in T], lower_bound=0.0)

# LINE FLOW CONSERVATION CONSTRAINTS:
# lines super source node sends 1 unit for each line
@constraint(
    master,
    lineBal[t in T],
    sum(d[params.nNodes+1, j, t] for j in lineDict[params.nNodes+1]) == 1
)

# line flows pass through stop nodes
@constraint(
    master,
    lineBal2[i in N, t in T],
    +sum(d[i, j, t] for j in lineDict[i]) -
    sum(d[j, i, t] for j in destination[i]) - d[params.nNodes+1, i, t] == 0
)

# line super sink node's flow demand is satisfied
@constraint(
    master,
    lineBal3[t in T],
    -sum(d[i, params.nNodes+2, t] for i in N) == -1
)

# line directions are already works bidirectional
@constraint(
    master,
    bidirect[i in N, j in destination[i], t in T],
    d[i, j, t] + d[j, i, t] <= 1
)

# only one unit of line flow outs from a node in a line
@constraint(
    master,
    outFlow[i in N, t in T],
    sum(d[i, j, t] for j in destination[i]) <= 1
)

# only one unit of line flow outs from a node in a line
@constraint(
    master,
    outFlow2[j in N, t in T],
    sum(d[i, j, t] for i in destination[j]) <= 1
)

# limits the number of stops in a line
@constraint(
    master,
    limitItinerary[t in T],
    sum(d[i, j, t] for i in N, j in destination[i]) == params.itinerary - 1
)

# subtour elimination by MTZ
# MTZs overstrains MASTER, no need
@constraint(master, tour[i in N, j in destination[i], t in T],
                            u[i,t] - u[j,t] + (params.itinerary-1)*(d[i,j,t]) <= params.itinerary-2)

@constraint(master, tour2[i in N, t in T],
                            u[i,t] <= params.itinerary-1)
#
# provide feasibility of the subproblem in the master problem
@constraint(master, feasibility[i in N],
                            sum(d[i,j,t] for j in destination[i], t in T)
                          + sum(d[j,i,t] for j in destination[i], t in T)
                            >= 1)
#
@objective(master, Min, p)

#=
fix(d[1,2,1], 1)
fix(d[2,5,1], 1)
fix(d[5,4,1], 1)
fix(d[4,6,1], 1)
fix(d[6,8,1], 1)
fix(d[8,10,1], 1)
fix(d[10,11,1], 1)

fix(d[9,15,2], 1)
fix(d[15,8,2], 1)
fix(d[8,10,2], 1)
fix(d[10,14,2], 1)
fix(d[14,13,2], 1)
fix(d[13,11,2], 1)
fix(d[11,12,2], 1)

fix(d[10,7,3], 1)
fix(d[7,15,3], 1)
fix(d[15,6,3], 1)
fix(d[6,3,3], 1)
fix(d[3,2,3], 1)
fix(d[2,4,3], 1)
fix(d[4,12,3], 1)

fix(d[12,11,4], 1)
fix(d[11,10,4], 1)
fix(d[10,8,4], 1)
fix(d[8,6,4], 1)
fix(d[6,3,4], 1)
fix(d[3,2,4], 1)
fix(d[2,1,4], 1)
=#

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
    set_optimizer_attribute(subProb,"Presolve", 0)
    set_optimizer_attribute(subProb,"Heuristics",0.0)
    set_optimizer_attribute(subProb,"Cuts",0)
    #set_optimizer_attribute(subProb,"Method",0)

    # dual variable for primal problem pTran's passenger flow conservation constraints - free variable
    @variable(subProb, lambda[i in N, k in passenger[i]])

    # dual variable for primal problem pTran's linking constraints with d_ijt
    @variable(subProb, mu[i in N, j in destination[i], t in T], lower_bound=0)

    # dual variable for primal problem pTran's transfer constraints
    @variable(subProb, teta[i in N, j in destination[i], k in passenger[i], t in T], lower_bound=0)

    # dual variable for primal problem pTran's maximum itinerary constraint - free variable
    @variable(subProb, beta[t in T])

    @variable(subProb, dual_artificial[i in N, j in destination[i], t in T])

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
                                + beta[t] + dual_artificial[i,j,t] == 0)

    # objective term for  primal problem pTran's passenger flow conservation constraint
    ex1 = @expression(subProb, dualObj1, sum(rhs[i,k]*lambda[i,k] for i in N, k in passenger[i]))

    # objective term for primal problem pTran's route itinerary constraint
    ex2 = @expression(subProb, dualObj2, sum((params.itinerary-1)*beta[t] for t in T))

    # objective term for primal Problme pTran's capacity-link constraints
    ex3 = @expression(subProb, dualObj3, sum(dual_artificial[i,j,t]*d_current[i,j,t]
                                              for i in N, j in destination[i], t in T))

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

        #construct terms in the optimality cut
        ex4 = @expression(master, sum(rhs[i,k]*value(lambda[i,k]) for i in N, k in passenger[i]))
        ex5 = @expression(master, sum((params.itinerary-1)*value(beta[t]) for t in T))
        ex6 = @expression(master, sum(value(dual_artificial[i,j,t])*d[i,j,t] for i in N, j in destination[i], t in T))

        @constraint(master, p >= ex4 + ex5 + ex6 )
    end

    if t_status_sub == MOI.DUAL_INFEASIBLE
        println(
            "\nThere is an  extreme ray, adding the corresponding constraint",
        )

        # in the backend we will retrieve the extreme rays of unbounded subproblem cases
        grb = backend(subProb)

        # to store extreme rays corresponding the subproblem
        lambdaF = zeros((params.nNodes, params.nNodes))
        muF = zeros((params.nNodes, params.nNodes, params.nLine))
        tetaF = zeros((params.nNodes, params.nNodes, params.nNodes, params.nLine))
        betaF = zeros((params.nLine))

        for i in N, j in destination[i], k in passenger[i], t in T
            r = Ref{Cdouble}()
            a = lambda[i,k].index.value
            Gurobi.GRBgetdblattrelement(grb, GRB_DBL_ATTR_UNBDRAY,a,r)
            lambdaF[i,k] = r.x

            r = Ref{Cdouble}()
            a = beta[t].index.value
            Gurobi.GRBgetdblattrelement(grb, GRB_DBL_ATTR_UNBDRAY,a,r)
            betaF[t] = r.x

            r = Ref{Cdouble}()
            a = dual_artificial[i,j,t].index.value
            Gurobi.GRBgetdblattrelement(grb, GRB_DBL_ATTR_UNBDRAY,a,r)
            dual_artificialF[i,j,t] = r.x
        end

        ex7 = @expression(master, sum(rhs[i,k]*lambdaF[i,k] for i in N, k in passenger[i]))
        ex8 = @expression(master, sum((params.itinerary-1)*betaF[t] for t in T))
        ex9 = @expression(master, sum(dual_artificialF[i,j,t]*d[i,j,t]
                                  for i in N, j in destination[i], t in T))

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
