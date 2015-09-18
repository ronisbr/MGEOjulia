#==#
#
# @brief Check the dominance of a Pareto point.
# @author Ronan Arraes Jardim Chagas
# @date 2015-04-03
#
# @param[in] mgeoData Structure with the configuration of MGEO.
# @param[in] candidatePoint Candidate point to be added in the frontier.
# @param[inout] paretoFrontier Pareto frontier.
#
# @retval TRUE The candidate point was added to the list.
# @retval FALSE The candidate point was not added to the list.
#
#==#

function checkDominance(mgeoData::MGEOStructure,
                        candidatePoint::ParetoPoint,
                        paretoFrontier::Array{ParetoPoint,1},
                        )
    addPoint = false

    # If Pareto frontier is empty, then return true.
    if (length(paretoFrontier) == 0)
        push!(paretoFrontier, candidatePoint)
        return true
    end

    remove = Int64[]

    # Loop through the entire list of Pareto points.
    for p=1:length(paretoFrontier)
        # Variable to check if the point in the Pareto frontier list is
        # dominated by the candidate.
        pointDominated = true
        candidateDominated = true

        # Loop through objective functions.
        for i=1:mgeoData.nf
            # If both objective functions are equal (given mgeoEps), then
            # nothing can be stated about the dominance of both points.
            if( abs(candidatePoint.f[i] - paretoFrontier[p].f[i]) <
               mgeoData.mgeoEps )
                continue
            end

            if( candidatePoint.f[i] < paretoFrontier[p].f[i] )
                candidateDominated = false
            elseif ( candidatePoint.f[i] > paretoFrontier[p].f[i] )
                pointDominated = false
            end
        end

        # If the point is dominated by the candidate and dominates the candidate
        # at the same time, then the candidate and the point is the same. Thus,
        # mark to do not add the candidate and exit the loop.

        if (pointDominated && candidateDominated)
            addPoint = false
            break
        end

        # If the point in the Pareto frontier is dominated, then mark it to be
        # excluded.
        if (pointDominated)
            push!(remove, p)
        end

        # If the candidate is dominated by any point in the list, stop the
        # search and do not add it.
        if (candidateDominated)
            addPoint = false
            break
        # Otherwise, add the point to the list.
        else
            addPoint = true
        end
    end

    # Remove the points dominated by the candidate point.
    for i=1:length(remove)
        deleteat!(paretoFrontier, remove[i]-(i-1))
    end

    # Check if the point must be added.
    if(addPoint)
        push!(paretoFrontier, candidatePoint)
    end

    addPoint
end

#==#
#
# @brief Configure the design variables.
# @author Ronan Arraes Jardim Chagas
# @date 2015-04-15
#
# @param[in] n The number of design variables.
# @param[in] bits The number of bits for each variable.
# @param[in] min The minimum value for all design variables.
# @param[in] max The maximum value for all design variables.
#
#==#

function confDesignVars(n::Int64,
                        bits::Int64,
                        min::Float64,
                        max::Float64)
    # Check if the number of bits are larger than 0.
    ( bits <= 0 ) && throw(MGEOArgumentError("The number of bits must not be 0."))

    # Check if the min value is smaller than max.
    ( min >= max ) && throw(MGEOArgumentError("The minimum value must be less than the maximum value for a variable."))

    # Create the array of design variables.
    designVariables = Array(DesignVariable, n)

    # Full scale for each variable.
    fullscale = (1 << bits) - 1

    # Factors to convert the string to real numbers.
    factors = 2.^collect(0:1:(bits-1))

    for i=1:n
        designVariables[i] = DesignVariable(bits,
                                            min,
                                            max,
                                            factors,
                                            fullscale,
                                            (i-1)*bits+1,
                                            "Var. $i")
    end

    designVariables
end

#==#
#
# @brief Configure the design variables.
# @author Ronan Arraes Jardim Chagas
# @date 2015-04-15
#
# @param[in] bits List containing the number of bits for each design variable.
# @param[in] min List containing the minimum of each design variable.
# @param[in] max List containing the maximum of each design variable.
# @param[in] varNames List containing the name of each design variable.
#
#==#

function confDesignVars{S<:String}(bits::Array{Int64,1},
                                   min::Array{Float64,1},
                                   max::Array{Float64,1},
                                   varNames::Array{S,1})
    # Check if the size of arrays is correct.
    n = length(bits)

    ( ( length(min) != n ) ||
     ( length(max) != n ) ||
     ( length(varNames) != n ) ) && throw(ArgumentError)

    # Number of bits.
    numBits = 0

    # Create the array of design variables.
    designVariables = Array(DesignVariable, n)

    for i = 1:n
        # Check if the number of bits are larger than 0.
        ( bits[i] <= 0 ) && throw(MGEOArgumentError("The number of bits must not be 0."))

        # Check if the min value is smaller than max.
        ( min[i] >= max[i] ) && throw(MGEOArgumentError("The minimum value must be less than the maximum value for a variable."))

        designVariables[i] = DesignVariable(bits[i],
                                            min[i],
                                            max[i],
                                            2.^collect(0:1:(bits[i]-1)),
                                            (uint64(1) << bits[i]) - 1,
                                            numBits+1,
                                            varNames[i])
        numBits += bits[i]
    end

     designVariables
end

#==#
#
# @brief Configure the design variables.
# @author Ronan Arraes Jardim Chagas
# @date 2015-04-15
#
# @param[in] bits List containing the number of bits for each design variable.
# @param[in] min List containing the minimum of each design variable.
# @param[in] max List containing the maximum of each design variable.
#
#==#

function confDesignVars(bits::Array{Int64,1},
                        min::Array{Float64,1},
                        max::Array{Float64,1})
    # Create an array with variable names.
    varNames = Array(String, length(bits))

    for i=1:length(varNames)
        varNames[i] = "Var. $i"
    end

    confDesignVars(bits, min, max, varNames)
end

#==#
#
# @brief Convert the string to real values.
# @author Ronan Arraes Jardim Chagas
# @date 2015-04-15
#
# @param[in] string String.
# @param[in] designVars Array with the description of design variables.
#
# @return Array(Float64,n) with the real value of the design variables, where n
# is the number of design variables.
#
#==#

function convertStringToNumber(designVars::Array{DesignVariable,1},
                               string::BitArray{1})
    # Number of variables.
    n = length(designVars)

    # Array of variables.
    vars = Array(Float64, n)
    # Loop for each design variable.
    for i = 1:n
        # Convert the bits to integer.
        string_index_i = designVars[i].index
        string_index_f = designVars[i].index+designVars[i].bits-1
        varInt =
            designVars[i].factors'*string[string_index_i:string_index_f]

        # Convert the integer to real number.
        vars[i] = designVars[i].min + (designVars[i].max - designVars[i].min)*
            float(varInt[1])/float(designVars[i].fullScale)
    end

    vars
end

#==#
#
# @brief Convert the structure that stores the Pareto frontier to an array.
# @author Ronan Arraes Jardim Chagas
# @date 2015-04-30
#
# @param[in] paretoFrontier Pareto frontier.
#
# @retun Array{Float64,2} in which the i-th line of the array is
#     [var[1] var[2] ... var[n] f[1] f[2] ... f[nf]],
# where n is the number of design variables and nf is the number of objective
# functions.
#
#==#

function pfToArray(paretoFrontier::Array{ParetoPoint,1})
    # Get the number of points.
    np = length(paretoFrontier)

    # Get the number of variables.
    n = length(paretoFrontier[1].vars)

    # Get the number of objective functions.
    nf = length(paretoFrontier[1].f)

    # Allocate the array.
    pfArray = zeros(np, n+nf)

    # Fill the array.
    for i=1:np
        pfArray[i,:] = [paretoFrontier[i].vars' paretoFrontier[i].f']
    end

    pfArray
end

#==#
#
# @brief Save the Pareto frontier to a CSV file.
# @author Ronan Arraes Jardim Chagas
# @date 2015-05-03
#
# @param[in] filename File name.
# @param[in] paretoFrontier Pareto frontier.
# @param[in] designVars Array with the description of design variables.
#
# @remarks The output file will have the following structure
#
#     var[1],var[2],...,var[n],f[1],f[2],...,f[nf],
#
# where n is  the number of design variables and nf is the number of objective
# functions.
#
#==#

function savePFtoCSV{S<:String}(filename::S,
                                designVars::Array{DesignVariable,1},
                                paretoFrontier::Array{ParetoPoint,1})
    # Convert Pareto frontier.
    pfArray = pfToArray(paretoFrontier)

    # Number of variable names.
    N = length(designVars)

    # Number of objective functions.
    nf = size(pfArray,2) - N

    # Array of variable names.
    varNames = Array(String, N)

    for i = 1:N
        varNames[i] = designVars[i].name
    end

    # Array of function names.
    funcNames = Array(String, nf)

    for i = 1:nf
        funcNames[i] = "Obj. Func. $i"
    end

    # Write the information to the file.
    writedlm(filename, [ [varNames' funcNames']; pfArray], ",")
end

#==#
#
# @brief Sort the points in the Pareto Frontier.
# @author Ronan Arraes Jardim Chagas
# @date 2015-04-30
#
# @param[inout] paretoFrontier Pareto frontier.
# @param[in] fobj Function to sort.
#
#==#

function sortParetoFrontier(paretoFrontier::Array{ParetoPoint,1}, fobj::Int64)
    # Check the parameter.
    nf = length(paretoFrontier[1].f)

    ( (fobj < 1) || (fobj > nf) ) && throw(MGEOArgumentError("fobj must be between 1 and the number of objective functions."))

    sort!(paretoFrontier, lt=(a,b)->begin
        # Check if data is valid.
        return (a.f[fobj] < b.f[fobj])
    end
          )
end
