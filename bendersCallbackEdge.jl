"""
MASTER IS edge formulation

BENDERS WITH CALLBACK
"""

flowq = ODistance(true,false)
pth = ODistance(false,true)
sup, rhs = rightHandSide()
params = problemParameters()
N, NT, T, destination, passenger,lineDict = requiredDictionaries()

#params.bigNumber = 15570
# subproblem model
function build_SubProb()

  # we use direct_model function since we will deal with the backend of the model
  # in the backend we will retrieve the extreme rays of unbounded subproblem cases
  subProb = direct_model(Gurobi.Optimizer())

    # we need info about the unbounded subproblem
    set_optimizer_attribute(subProb,"InfUnbdInfo",1)
    set_optimizer_attribute(subProb,"Presolve",2)
    #set_optimizer_attribute(subProb,"Heuristics",0.0)
    #set_optimizer_attribute(subProb,"Cuts",0)
    #set_optimizer_attribute(subProb,"Method",2)
    #set_optimizer_attribute(subProb,"FeasibilityTol", 1e-2)

    # dual variable for primal problem pTran's passenger flow conservation constraints - free variable
    @variable(subProb, lambda[i in N, k in passenger[i]])

    # dual variable for primal problem pTran's linking constraints with d_ijt
    @variable(subProb, mu[i in N, j in destination[i], t in T; i<j], upper_bound=0)

    # dual variable for primal problem pTran's transfer constraints
    @variable(subProb, teta[i in N, j in destination[i], k in passenger[i], t in T; i!=k], upper_bound=0)

    # dual artificial variable for the third term of the subproblems objective function
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
                                 + artificial[i,j,t] == 0)

    # constraint for primal problem pTran's transfer variable
    @constraint(subProb, dualTransfer[i in N, j in destination[i], k in passenger[i], t in T; i!=k],
                                - teta[i,j,k,t] <= params.penalty)


    # objective term for  primal problem pTran's passenger flow conservation constraint
    ex1 = @expression(subProb, dualObj1, sum(rhs[i,k]*lambda[i,k] for i in N, k in passenger[i]))

  # dualObj3 can be added if it is needed
  return (subProb, lambda, artificial, ex1)

end

master = direct_model(Gurobi.Optimizer())

(subProb, lambda, artificial, ex1) = build_SubProb()

set_optimizer_attribute(master,"Presolve",0)
set_optimizer_attribute(master,"Heuristics",0.0)
#set_optimizer_attribute(master,"Cuts",3)
#set_optimizer_attribute(master, "Threads",1) # sometimes parallel cores find solutions in other nodes so lead to wrong answers
#set_optimizer_attribute(master,"Method",0)

# create a surrogate variable to hold subProblem's objective in the master
@variable(master, p, lower_bound = -10e6)

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

# place specified number of stops in a route
@constraint(master, limitItinerary[t in T],
                        sum(d[i,j,t] for i in N, j in destination[i] if i<j) == params.itinerary-1)

# place each stop at least one time in the route plan
@constraint(master, atleastOne[i in N], sum(z[i,t] for t in T) >= 1)

@objective(master, Min, p)

#iter_num = 0
cb_calls = Cint[]

function benders(cb_data, cb_where::Cint)

  global iter_num
  #iter_num+=1

  #println("Iteration no : ", iter_num)

  # You can select where the callback is run
  if cb_where != GRB_CB_MIPSOL && cb_where != GRB_CB_MIPNODE
    return
  end

  #
  if cb_where == GRB_CB_MIPNODE
    resultP = Ref{Cint}()
    GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_STATUS, resultP)
    if resultP[] != GRB_OPTIMAL
      return  # Solution is something other than optimal.
    end
  end
  #

  # You can reference variables outside the function as normal
  push!(cb_calls, cb_where)

  # Before querying `callback_value`, you must call:
  Gurobi.load_callback_variable_primal(cb_data,cb_where)

  global d_current = callback_value.(cb_data, d)
  global master_current = callback_value(cb_data, p)

  # objective term for  primal problem master's linking constraints
  ex3 = @expression(subProb, sum(artificial[i,j,t] * d_current[i,j,t] for i in N, j in destination[i], t in T if i<j))

  @objective(subProb, Max, ex1 + ex3)

  # we will add the commodity flows to eliminate subtours from the master's soluton

  optimize!(subProb)

  #println("\n objective if subproblem", objective_function(subProb))

  global p_status = termination_status(subProb)
  global prim_status = primal_status(subProb)

  global sub_current = objective_value(subProb)

  println("\n master_CURRENT: ", master_current)




  # subtour detection and elimination
  if cb_where == GRB_CB_MIPSOL

    for t in T
      global g = SimpleGraph(params.nNodes+2)
      for i in NT, j in lineDict[i]
        if (i<j || i==params.nNodes+1) && d_current[i,j,t] > 0.5
          add_edge!(g, i,j)
        end
      end

      cyc = connected_components(g) # find connected segments in a route
      global cycles = Vector{Vector{Int64}}() # create for storing cycles in a route
      # seperate cycles from the graph by searchin conected components
      for i in 1:size(cyc,1)
        if length(cyc[i]) > 2 && params.nNodes+1 ∉ cyc[i]
          push!(cycles, cyc[i])
        end
      end
      for r in 1:length(cycles)
        subtour_cons = @build_constraint(sum(d[i,j,t]*d_current[i,j,t] for i in cycles[r], j in cycles[r], t in t if i<j && j in destination[i]) <= length(cycles[r])-1)
        MOI.submit(master, MOI.LazyConstraint(cb_data), subtour_cons)
        println("subtour elimination constraint: ", subtour_cons)
      end
    end
  end
  #





  if p_status == MOI.INFEASIBLE_OR_UNBOUNDED || p_status == MOI.DUAL_INFEASIBLE
    println("\nAdding a feasibility cut")

          # in the backend we will retrieve the extreme rays of unbounded subproblem cases
          grb = backend(subProb)

          # to store extreme rays corresponding the subproblem
          global lambdaF = zeros((params.nNodes,params.nNodes))
          global artificialF = zeros((params.nNodes,params.nNodes,params.nLine))

          for i in N, j in destination[i], k in passenger[i], t in T
            if i<j
              r = Ref{Cdouble}()
              a = lambda[i,k].index.value
              Gurobi.GRBgetdblattrelement(grb, GRB_DBL_ATTR_UNBDRAY,a,r)
              lambdaF[i,k] = r.x

              r = Ref{Cdouble}()
              a = artificial[i,j,t].index.value
              Gurobi.GRBgetdblattrelement(grb, GRB_DBL_ATTR_UNBDRAY,a,r)
              artificialF[i,j,t] = r.x
            end
          end

          ex7 = @expression(master, sum(rhs[i,k]*lambdaF[i,k] for i in N, k in passenger[i]))
          ex9 = @expression(master, sum(artificialF[i,j,t]*d[i,j,t]
                                                    for i in N, j in destination[i], t in T if i<j))


    feas_cons = @build_constraint(ex7 + ex9 <= 0)

    MOI.submit(master, MOI.LazyConstraint(cb_data), feas_cons)
    println("feasibility cut is added")
  end


  if prim_status == MOI.FEASIBLE_POINT && sub_current != master_current
    println("\nAdding an optimality cut")

    #construct terms in the optimality cut
    ex4 = @expression(master, sum(rhs[i,k]*value(lambda[i,k]) for i in N, k in passenger[i]))
    ex6 = @expression(master, sum(value(artificial[i,j,t]) * d[i,j,t] for i in N, j in destination[i], t in T if i<j))

    # create the optimality cut
    optim_cons = @build_constraint( p >= ex4 + ex6)

    # add the optimality cut
    MOI.submit(master, MOI.LazyConstraint(cb_data), optim_cons)
  end
  #

  if prim_status == MOI.FEASIBLE_POINT && abs(sub_current - master_current)<0.00001
    @info("SUBPROBLEM AND MASTER PROBLEM EQUAL")
    return
  end #

end

# You _must_ set this parameter if using lazy constraints.
MOI.set(master, MOI.RawParameter("LazyConstraints"), 1)
MOI.set(master, Gurobi.CallbackFunction(), benders)

#MOI.set(master,MOI.LazyConstraintCallback(), benders)

optimize!(master)
