
"""
load data from CSV files and determine the parameters of the problem

load the data instance here

"""


# import CSVs of OD and Distance matrices
function ODistance(flow::Bool, distance::Bool)
    # import distance and path data
    # file paths must be updated for different machines
    demand = CSV.read("/Users/akrm/Dropbox/Julia/transitdata/mandl/MandlDemand.csv", header = false, delim = ';', decimal=',', DataFrame)
    dist = CSV.read("/Users/akrm/Dropbox/Julia/transitdata/mandl/MandlTravelTimes.csv", header = false, delim = ';', decimal=',', DataFrame)

    #demand = CSV.read("/Users/akrm/Dropbox/Julia/transitdata/kayseri/FlowMatrix79Use.csv", header = false, delim = ';', decimal=',', DataFrame)
    #dist = CSV.read("/Users/akrm/Dropbox/Julia/transitdata/kayseri/Distance79Use.csv", header = false, delim = ';', decimal=',', DataFrame)


    # converting the dataframes into arrays and make some data preprocess
    flowq = Matrix(demand)
    pth = Matrix(dist)

    # round the float numbers
    flowq = round.(Int, flowq)
    pth = round.(Int, pth)



                # make the OD is one directional
                for i in 1:size(flowq,1),  j in 1:size(flowq,1)
                    if j>i
                        flowq[j,i] = flowq[j,i] + flowq[i,j]
                        flowq[i,j] = 0
                    end
                end
                #

    # distances are symmetric by diagonal
    for i = 1:size(pth,1), j=1:size(pth,1)
        pth[j,i] = pth[i,j]
    end
    if flow == true
        return flowq
    elseif distance == true
        return pth
    end
end

@enum Subtour MTZ MULTICOMMODITY PASSENGERFLOW

"""
set of parameters of the problem
"""
mutable struct Parameters
    "number of stops in the network"
    nNodes::Int
    "number of desired lines"
    nLine::Int
    "number of seats in the transport vehicle"
    cap::Vector{Int64}
    "transfer penalty"
    penalty::Float64
    "the big number of capacity constraint"
    bigNumber::Float64
    "number of nodes allowed in a line"
    itinerary::Int
    "subtour type"
    subtour::Subtour
    Parameters() = new()
end

function problemParameters()
    params = Parameters()
    params.nNodes = size(flowq,1)
    params.nLine = 4
    params.cap = fill(50,params.nLine)
    params.penalty = 5
    #params.bigNumber = round(551720/50)
    params.bigNumber = round(15570)
    params.itinerary = 8
    params.subtour = MTZ
    return params
end

"construct the right hand side of the passenger flow conversation"
function rightHandSide()
    params = problemParameters()
    flowq = ODistance(true,false)
    # create total supply for source nodes (1st layer)
    sup = sum(flowq[1:params.nNodes, 1:params.nNodes], dims=2)

    # right hand side values of passenger balance constraints
    rhs = zeros(Int32, params.nNodes, params.nNodes)
    for i in 1:params.nNodes, k in 1:params.nNodes
        rhs[i,k] = -flowq[k,i]
        rhs[k,k] = sup[k]
    end
    return sup, rhs
end

"create dictionaries for range of decision variable sets"
function requiredDictionaries()
    params = problemParameters()
    # define the sets of the problem
    N =  collect(1:params.nNodes)            # stop network
    NT = collect(1:params.nNodes+2)          # whole network
    T =  collect(1:params.nLine)             # lines

    # destination dictionary for "x" passenger flow variables
    destination = Dict(i => [] for i in N)
    for i in N
        for j in N
            if pth[i,j] > 0
                push!(destination[i], j)
            end
        end
    end

    # passenger types are assigned to stops
    passenger = Dict()
    for i in 1:length(N)
        passenger[i] = N
    end

    # destination matrix for line "d" line flow variables
    lineDict = Dict(i => [] for i in NT)
    for i in N
        for j in N
            if pth[i,j] > 0
                push!(lineDict[i], j)
            end
        end
        push!(lineDict[i], params.nNodes+2)         # adding line sink
    end
    for j in N
        push!(lineDict[params.nNodes+1], j)         # adding line source
    end

    return N, NT, T, destination, passenger, lineDict
end

# if params.subtour == MULTICOMMODITY
"add lines as nodes in the network for MULTICOMMODITY subtour"
function addLineNode()
            forbiddenDistance = 200
            # adding lines as nodes
            # arrange distance matrix according to the lines as nodes
            tem = fill(forbiddenDistance, nLine, nNodes)
            tem2 = transpose(tem)
            tem3 = zeros(nLine,nLine)
            tem4 = vcat(tem2, tem3)
            tem5 = vcat(pth, tem)
            pthh = hcat(tem5, tem4)
            pth2 = zeros(nNodes+2*nLine, nNodes+2*nLine, nLine)
            for i=1:nNodes+2*nLine, j=1:nNodes+2*nLine, t=1:nLine
                if i>nNodes && j<=nNodes && i == t+nNodes
                    pth2[i,j,t] = forbiddenDistance
                elseif i<=nNodes && j > nNodes && j==(t+nNodes)+nLine
                    pth2[i,j,t] = forbiddenDistance
                elseif i<=nNodes && j<=nNodes
                    pth2[i,j,t] = pthh[i,j]
                end
            end
        return pth2
end

# if the subtour elimination is handled with a MULTICOMMODITY flow variable
# add an extra dictionary for that variable
"create dictionary for range of MULTICOMMODITY subtour flow variable"
function additionalDictionaries()
        # destination matrix for 'y' variables
        lineDict2 = Dict()
        for t in T
            destination2 = Dict(i => [] for i in NT)
            for i in 1:size(pth2,1)
                for j in 1:size(pth2,1)
                    if pth2[i,j,t] > 0
                        push!(destination2[i], j)
                    elseif i == nNodes+nLine+t && j <=nNodes
                        push!(destination2[i], j)
                    end
                end
            end
            lineDict2[t] = destination2
        end

        return lineDict2
end
