"""
designed as a module file

the main stream
"""

using DataFrames, CSV, DelimitedFiles, LightGraphs
using JuMP, Gurobi, CPLEX

include("loadData.jl")

# arc formulation with subtour elimination MTZ
#include("arcMTZ.jl")

# edge formulation wothout subtour elimination constraints
#include("EdgeNoSubtour.jl")

# edge formulation with commodity subtour elimination constraints
#include("EdgeWithCommoditySubtour.jl")

#optimize!(pTran)


#include("exportResults.jl")
