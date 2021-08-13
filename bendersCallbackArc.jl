"""
MASTER IS ARC MTZ
BENDERS WITH CALLBACK
"""

flowq = ODistance(true,false)
pth = ODistance(false,true)
sup, rhs = rightHandSide()
params = problemParameters()
N, NT, T, destination, passenger,lineDict = requiredDictionaries()
# the new commodity set including line commodities
# K = [N; collect(params.nNodes+3:params.nNodes+3+params.nLine-1)]

#params.bigNumber = 15570
# subproblem model
function build_SubProb()

  # we use direct_model function since we will deal with the backend of the model
  # in the backend we will retrieve the extreme rays of unbounded subproblem cases
  subProb = direct_model(Gurobi.Optimizer())
  set_silent(subProb)

  # we need info about the unbounded subproblem
  set_optimizer_attribute(subProb,"InfUnbdInfo",1)
  set_optimizer_attribute(subProb,"Presolve",2)
  set_optimizer_attribute(subProb,"BarConvTol",0)
  set_optimizer_attribute(subProb,"ScaleFlag",2)
  set_optimizer_attribute(subProb,"DegenMoves",0)



  #set_optimizer_attribute(subProb,"Heuristics",0.0)
  #set_optimizer_attribute(subProb,"Cuts",0)
  set_optimizer_attribute(subProb,"Method",2)
  #set_optimizer_attribute(subProb,"FeasibilityTol", 1e-2)


  # dual variable for primal problem pTran's passenger flow conservation constraints - free variable
  @variable(subProb, lambda[i in N, k in passenger[i]])

  # dual variable for primal problem pTran's linking constraints with d_ijt
  @variable(subProb, mu[i in N, j in destination[i], t in T], upper_bound=0)

  # dual variable for primal problem pTran's transfer constraints
  @variable(subProb, teta[i in N, j in destination[i], k in passenger[i], t in T], upper_bound=0)

  # dual variable for primal problem pTran's maximum itinerary constraint - free variable
  #@variable(subProb, beta[t in T])

  @variable(subProb, dual_artificial[i in N, j in destination[i], t in T])

  # constraint for primal problem pTran's passenger flow variable
  @constraint(subProb, dualFlow[i in N, j in destination[i], k in passenger[i], t in T; i==k],
                              lambda[i,k] - lambda[j,k]
                              - sum(teta[j,r,k,t] for r in destination[j] if r!=i)
                              + mu[i,j,t] <= pth[i,j])

  # constraint for primal problem pTran's passenger flow variable
  @constraint(subProb, dualFlow2[i in N, j in destination[i], k in passenger[i], t in T; i!=k],
                              lambda[i,k] - lambda[j,k]
                              + teta[i,j,k,t]
                              - sum(teta[j,r,k,t] for r in destination[j] if r!=i)
                              + mu[i,j,t] <= pth[i,j])

  # constraint for primal problem pTran's transfer variable
  @constraint(subProb, dualTransfer[i in N, j in destination[i], k in passenger[i], t in T],
                              - teta[i,j,k,t] <= 5)

  # constraint for primal problem pTran's line flow variable (d[i,j,t])
  @constraint(subProb, dualLink[i in N, j in destination[i], t in T],
                              - params.bigNumber*mu[i,j,t] - params.bigNumber*mu[j,i,t]
                              #+ beta[t]
                              + dual_artificial[i,j,t] == 0)

  # objective term for  primal problem pTran's passenger flow conservation constraint
  ex1 = @expression(subProb, dualObj1, sum(rhs[i,k]*lambda[i,k] for i in N, k in passenger[i]))

  # dualObj3 can be added if it is needed
  return (subProb, lambda, mu, teta, dual_artificial, ex1)

end

master = direct_model(Gurobi.Optimizer())

(subProb, lambda, mu, teta, #=beta,=# dual_artificial, ex1#=, ex2=# ) = build_SubProb()

set_optimizer_attribute(master,"Presolve",0)
set_optimizer_attribute(master,"Heuristics",0.0)
#set_optimizer_attribute(master,"Cuts",3)
#set_optimizer_attribute(master, "Threads",1) # sometimes parallel cores find solutions in other nodes so lead to wrong answers
set_optimizer_attribute(master,"Method",0)

# create a surrogate variable to hold subProblem's objective in the master
@variable(master, p, lower_bound = -1e6)

# line flow decision variable
@variable(master, d[i in NT, j in lineDict[i], t in T], Bin)

# adding MTZ subtour elimination variable
@variable(master, u[i in N, t in T], lower_bound = 0.0)

# LINE FLOW CONSERVATION CONSTRAINTS:
#lines super source node sends 1 unit for each line
@constraint(master, lineBal[t in T],
                            sum(d[params.nNodes+1,j,t] for j in lineDict[params.nNodes+1]) == 1)

# line flows pass through stop nodes
@constraint(master, lineBal2[i in N, t in T],
                            + sum(d[i,j,t] for j in lineDict[i])
                            - sum(d[j,i,t] for j in destination[i])
                            - d[params.nNodes+1, i, t] == 0)

# line super sink node's flow demand is satisfied
@constraint(master, lineBal3[t in T],
                            - sum(d[i,params.nNodes+2,t] for i in N) == -1)

# line directions are already works bidirectional
@constraint(master, bidirect[i in N, j in destination[i], t in T],
                            d[i,j,t] + d[j,i,t] <= 1)

# only one unit of line flow outs from a node in a line
@constraint(master, outFlow[i in N, t in T],
                            sum(d[i,j,t] for j in destination[i]) <= 1)

# only one unit of line flow outs from a node in a line
@constraint(master, outFlow2[j in N, t in T],
                            sum(d[i,j,t] for i in destination[j]) <= 1)

# limits the number of stops in a line
@constraint(master, limitItinerary[t in T],
                            sum(d[i,j,t] for i in N, j in destination[i]) == params.itinerary-1)

# subtour elimination by MTZ
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


#iter_num = 0
cb_calls = Cint[]

function benders(cb_data, cb_where::Cint)

  global iter_num
  #iter_num+=1

  #println("Iteration no : ", iter_num)

  # You can select where the callback is run
  if cb_where != GRB_CB_MIPSOL #  && cb_where != GRB_CB_MIPNODE
      return
  end

  #=
  if cb_where == GRB_CB_MIPNODE
    resultP = Ref{Cint}()
    GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_STATUS, resultP)
    if resultP[] != GRB_OPTIMAL
      return  # Solution is something other than optimal.
    end
  end
  =#

  # You can reference variables outside the function as normal
  push!(cb_calls, cb_where)

  # Before querying `callback_value`, you must call:
  Gurobi.load_callback_variable_primal(cb_data,cb_where)

  d_current = callback_value.(cb_data, d)
  #d_current = round.(d_current)
  master_current = callback_value(cb_data, p)

  #println(d_current)

  # objective term for  primal problem master's linking constraints
  ex3 = @expression(subProb, sum(dual_artificial[i,j,t] * d_current[i,j,t] for i in N, j in destination[i], t in T))

  @objective(subProb, Max, ex1 + ex3)

  optimize!(subProb)

  #println("\n objective if subproblem", objective_function(subProb))

  p_status = termination_status(subProb)
  prim_status = primal_status(subProb)

  sub_current = objective_value(subProb)

  #=
  if objective_value(subProb) < 163210 && prim_status == MOI.FEASIBLE_POINT
    open("errorrrr.txt", "w") do f
      println(f, d_current)
    end
  end
  =#
  println("SUB_STATUS", p_status)
  println("\n sub_current: ", sub_current)
  println("\n master_CURRENT:  ", master_current)

  if prim_status == MOI.FEASIBLE_POINT && sub_current == master_current
    @info("MASTER AND SUBPROBLEM IS SAME")
  end

  if prim_status == MOI.FEASIBLE_POINT && sub_current > master_current
    println("\nAdding an optimality cut")

    #construct terms in the optimality cut
    ex4 = @expression(master, sum(rhs[i,k]*value(lambda[i,k]) for i in N, k in passenger[i]))
    #ex5 = @expression(master, sum((params.itinerary-1)*value(beta[t]) for t in T))
    ex6 = @expression(master, sum(value(dual_artificial[i,j,t])*d[i,j,t] for i in N, j in destination[i], t in T))

    # create the optimality cut
    optim_cons = @build_constraint( p >= ex4 + ex6)

    # add the optimality cut
    MOI.submit(master, MOI.LazyConstraint(cb_data), optim_cons)

    #println(optim_cons)

  end
  #
  if p_status == MOI.INFEASIBLE_OR_UNBOUNDED || p_status == MOI.DUAL_INFEASIBLE
    println("\nAdding a feasibility cut")

    # in the backend we will retrieve the extreme rays of unbounded subproblem cases
    grb = backend(subProb)

    # to store extreme rays corresponding the subproblem
    lambdaF = zeros((params.nNodes,params.nNodes))
    betaF = zeros(params.nLine)
    dual_artificialF = zeros((params.nNodes,params.nNodes,params.nLine))

    for i in N, j in destination[i], k in passenger[i], t in T
      r = Ref{Cdouble}()
      a = lambda[i,k].index.value
      Gurobi.GRBgetdblattrelement(grb, GRB_DBL_ATTR_UNBDRAY,a,r)
      lambdaF[i,k] = r.x

      r = Ref{Cdouble}()
      a = dual_artificial[i,j,t].index.value
      Gurobi.GRBgetdblattrelement(grb, GRB_DBL_ATTR_UNBDRAY,a,r)
      dual_artificialF[i,j,t] = r.x
    end

    ex7 = @expression(master, sum(rhs[i,k]*lambdaF[i,k] for i in N, k in passenger[i]))
    ex9 = @expression(master, sum(dual_artificialF[i,j,t]*d[i,j,t]
                                  for i in N, j in destination[i], t in T))

    feas_cons = @build_constraint( ex7 + ex9 <= 0)


    MOI.submit(master, MOI.LazyConstraint(cb_data), feas_cons)

  end
end

# You _must_ set this parameter if using lazy constraints.
MOI.set(master, MOI.RawParameter("LazyConstraints"), 1)
MOI.set(master, Gurobi.CallbackFunction(), benders)

#MOI.set(master,MOI.LazyConstraintCallback(), benders)

optimize!(master)




#=
t_status = termination_status(master)
p_status = primal_status(master)

if p_status == MOI.INFEASIBLE_POINT
    println("The problem is infeasible")
end

if t_status == MOI.INFEASIBLE_OR_UNBOUNDED
    master_current = 1000
    d_current = ones(204)
end
=#

#=
if p_status == MOI.FEASIBLE_POINT
    master_current = value(p)
    d_current = value.(d)
end

#
println("Status of the master problem is ", t_status,
        "\nwith master_current = ", master_current)
#

#
for i in N, j in destination[i], t in T
  if d_current[i,j,t]>0
    println(i,",",j,",",t, "= ",d_current[i,j,t])
  end
end
=#

#=
for i in N, j in destination[i], k in passenger[i], t in T
  if value(mu[i,j,k,t])>0
    println(i,",",j,",",k,",",t, "= ",value(mu[i,j,k,t]))
  end
end
=#

#=
for i in N, k in passenger[i]
  if value(lambda[i,k])>0
    println(i,",",k, "=",value(lambda[i,k]))
  end
end
=#
#@test termination_status(master) == MOI.OPTIMAL
#@test primal_status(master) == MOI.FEASIBLE_POIN



## PRIMAL PROBLEM
#=
#
primal = Model(optimizer_with_attributes(Gurobi.Optimizer, "Presolve" => 0))


# passenger flow decision variable
@variable(primal, x[i in N, j in destination[i], k in passenger[i], t in T], lower_bound=0)

@variable(primal, d[i in N, j in destination[i], t in T])

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
@constraint(primal, capacity[i in N, j in destination[i], k in passenger[i], t in T],
                            x[i,j,k,t] <= params.bigNumber*d[i,j,t])

@constraint(primal, capacity2[i in N, j in destination[i], k in passenger[i], t in T],
                            x[j,i,k,t] <= params.bigNumber*d[i,j,t])

#

# passenger flow cost between two stops - objective expression
@expression(primal, costobj1[i in N, j in destination[i], k in passenger[i], t in T],
                            pth[i,j]*x[i,j,k,t])

# transfer cost
@expression(primal, costobj2[i in N, j in destination[i], k in passenger[i], t in T],
                            params.penalty*trnsfer[i,j,k,t])

@constraint(primal, link4[t in T],
                        sum(d[i,j,t] for i in N,  j in destination[i])
                        <= params.itinerary-1)


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


# objective function
@objective(primal, Min, sum(costobj1) + sum(costobj2))
=#

#=
using Dualization

dual_model = dualize(primal; dual_names = DualNames("dual", ""))

print((dual_model))

optimize!(primal)
=#
#
