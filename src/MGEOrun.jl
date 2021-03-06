"""
### function MGEOrun(mgeoData::MGEOStructure, f_obj::Function, showDebug::Bool = false)

Run the MGEO.

##### Args

* mgeoData: Structure with the configuration of MGEO.
* f_obj: Objective function.
* showDebug: Print debug information.

##### Remarks

The objective function must have the following signature

        function f_obj(vars)
            1) Compute f[1], f[2], ..., f[nf] using vars[1], vars[2], ... var[n]
            2) return (validPoint, f)
        end

where nf is the number of objective functions, n is the number of design
variables, and validPoint is a Boolean value that indicates if vars yield to a
valid point.

##### Info

* **Author**: Ronan Arraes Jardim Chagas
* **Date**: 2015-04-30

"""

function MGEOrun(mgeoData::MGEOStructure,
                 f_obj::Function,
                 showDebug::Bool = false)
    # Pareto frontier.
    paretoFrontier = ParetoPoint[]

    # Maximum number of generations per run.
    ngenMaxPerRun = floor(mgeoData.ngenMax/mgeoData.runMax)

    # Get the number of design variables.
    numDesignVars = length(mgeoData.designVars)

    # Get the number of bits the string must have by summing the bits for each
    # design variable.
    numBits = 0
    for i = 1:numDesignVars
        numBits += mgeoData.designVars[i].bits
    end

    # Create the string.
    string = BitArray(numBits)

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    #                   Initialization of Pareto Frontier
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Initialize the string with one valid point in the Pareto Frontier.
    for run = 1:mgeoData.runMax
        # Sample a new string.
        rand!(string)

        # Convert string to real numbers.
        vars = convertStringToNumber(mgeoData.designVars, string)

        # Call the objective functions.
        (valid, f) = f_obj(vars)

        # Check if the number of objective functions are valid.
        (length(f) != mgeoData.nf) && throw(MGEOArgumentError("The number of objective functions returned by f_obj is not equal to the one specified in mgeoData."))

        if (valid)
            # Add the point to the Pareto Frontier.
            push!(paretoFrontier, ParetoPoint(vars, f))
            break
        end
    end

    # Check if a valid point was found.
    if (length(paretoFrontier) == 0)
        return paretoFrontier
    end

    # String.
    string = BitArray(numBits)

    # Number of available processes.
    np = nprocs()
    i  = 1
    nextidx() = (idx = i; i+=1; idx)

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    #                       Loop - Independent Runs
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    for run = 1:mgeoData.runMax
        if (showDebug)
            println("--------------------------------------------------------------------------------")
            println("RUN = $run")
            println("")
            println("Generations ($ngenMaxPerRun):")
            print("    0%")
        end

        # Sample a new string if it is not the first run.
        (run > 1) && rand!(string)

        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        #                   Loop - MGEO Generations
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

        genDebug = 0

        for ngenPerRun = 1:ngenMaxPerRun
            if (showDebug)
                genDebug += 1

                if (genDebug >= 0.1*ngenMaxPerRun)
                    genDebug = 0
                    print("...$(ngenPerRun/ngenMaxPerRun*100)%")
                end
            end

            # Choose which objective function will be used to compute the
            # adaptability and to assemble the rank.
            chosenFunc = rand(1:mgeoData.nf)

            # List of all points created after flipping the bits.
            candidatePoints = Array(ParetoPoint,numBits)

            # Array to sort the candidate points.
            fRank = Array(sRank,numBits)

            i = 1

            # Evaluate the objective functions using parallel computing.
            @sync begin
                for p = 1:np
                    if p != myid() || np == 1
                        @async begin
                            while true
                                idx = nextidx()

                                if idx > numBits
                                    break
                                end

                                string[idx] = !string[idx]
                                vars = convertStringToNumber(
                                    mgeoData.designVars, string)

                                # Evaluate the objective functions.
                                (valid, f) = remotecall_fetch(f_obj, p, vars)

                                # Check if the solution is valid.
                                if (valid)
                                    # Create the candidate point.
                                    candidatePoint = ParetoPoint(vars, f)

                                    # Add the result to the rank.
                                    fRank[idx] = sRank(true, idx, f[chosenFunc])
                                    candidatePoints[idx] = candidatePoint
                                else
                                    fRank[idx] = sRank(false, idx, 0.0)
                                end
                            end
                        end
                    end
                end
            end

            # Add the points to the Pareto frontier.
            numValidPoints = 0

            for i=1:length(candidatePoints)
                if (fRank[i].valid)
                    numValidPoints += 1
                    checkDominance(mgeoData, candidatePoints[i], paretoFrontier)
                end
            end

            # If there are no valid points, an independent run must be called.
            if (numValidPoints == 0)
                break
            end

            # Ranking.
            # ------------------------------------------------------------------
            # The points will be ranked according to the selected objective
            # function (chosenFunc).
            #
            # Notice that the invalid points will be placed after the valid
            # ones.

            sort!(fRank, lt=(a,b)->begin
                # Check if data is valid.
                if (a.valid && b.valid)
                    return a.f < b.f
                elseif (a.valid)
                    return true
                else
                    return false
                end
            end
                  )

            # Choose a bit to be changed for the new generation.
            bitAccepted = false

            while (bitAccepted == false)
                # Sample a valid bit in the string.
                b_sample = rand(1:numValidPoints)

                # Accept the change with probability r^(-tal), where r is the
                # rank of the bit.
                Pk = b_sample^(-mgeoData.tau)

                if (rand() <= Pk)
                    # The bit is accepted, then exit the loop.
                    string[fRank[b_sample].index] =
                        !string[fRank[b_sample].index]
                    bitAccepted = true
                end
            end
        end

        if (showDebug)
            println("")
            println("")
            println("Pareto frontier information:")
            println("    Number of points: $(length(paretoFrontier))")
            println("")
        end
    end

    paretoFrontier
end
