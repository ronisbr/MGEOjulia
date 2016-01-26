module MGEOjulia

export DesignVariable, ParetoPoint, MGEOStructure
export confDesignVars, pfToArray, sortParetoFrontier, savePFtoCSV
export MGEOrun

################################################################################
#                                  Structures
################################################################################

"""
### type DesignVariable

Structure that defines the limits of the design variables.

"""

type DesignVariable
    # Number of bits.
    bits::Int64
    # Minimum allowed value for the design variable.
    min::Float64
    # Maximum allowed value for the design variable.
    max::Float64
    # Factors to convert the string to real numbers.
    factors::Array{Int64,1}
    # Full scale of the design variable.
    fullScale::Int64
    # Index in the string.
    index::Int64
    # Name of the variable.
    name::AbstractString
end

"""
### type ParetoPoint

Structure that defines a point in the Pareto frontier.

"""

type ParetoPoint
    # Design variables.
    vars::Array{Float64,1}
    # Objective functions.
    f::Array{Float64, 1}
end

"""
### type sRank

Structure used to sort the candidate points.

"""

type sRank
    # Is the value valid?
    valid::Bool
    # Index of the flipped bit.
    index::Int64
    # Objective function value after flipping the bit.
    f::Float64
end

"""
### type MGEOStructure

Structure to store the configuration of MGEO.

"""

type MGEOStructure
    # Number of objective functions.
    nf::Int64
    # Parameter to set the search determinism.
    tau::Float64
    # Maximum number of generations.
    ngenMax::Int64
    # Maximum number of independent runs (reinitializations).
    runMax::Int64
    # Structure that store the configuration for each design variable.
    designVars::Array{DesignVariable,1}
    # Epsilon to compare objective functions.
    mgeoEps::Float64
end

################################################################################
#                                  Exceptions
################################################################################

# Exception: Invalid arguments.
type MGEOArgumentError <: Exception
    msg::AbstractString
end
Base.showerror(io::IO, e::MGEOArgumentError) =
    print(io, e.msg)

################################################################################
#                                    Files
################################################################################

include("MGEOaux.jl")
include("MGEOrun.jl")

end # module
