using Distributions
using UnicodePlots
using Printf

"""
    function randomLengths(numAtoms, averageLenght, dispersity)

Randomly generate a list of chain lenghts with a sum of exaxtly numAtoms.

# Arguments
 - `numAtoms`
 - `averageLenght`
 - `dispersity`: see https://en.wikipedia.org/wiki/Dispersity
"""
function randomLengths(numAtoms, averageLenght, dispersity)
    shape = 1/(dispersity-1)
    scale = averageLenght*(dispersity-1)
    dist = Gamma(shape, scale)

    lengths = Int64[]
    atomCount = 0

    while atomCount < numAtoms
        length = ceil(rand(dist))
        atomCount += length
        push!(lengths, length)
    end

    lengths[end] -= atomCount - numAtoms

    return lengths
end

function printChainLengthStatistics(chainLengths)
    numChains = length(chainLengths)
    numAtoms = sum(chainLengths)
    Mn = numAtoms/numChains
    Mw = sum(chainLengths.^2)/numAtoms
    dispersity = Mw/Mn

    statisticsString = @sprintf "n=%d   Σ=%d   μ=%.1f   Ð=%.2f" numChains numAtoms Mn dispersity

    println()
    histogram(chainLengths, width=60, title="Chain Length Distribution", xlabel=statisticsString, border= :bold, margin=10) |> show
    println("\n")
end

function printAngleStatistics(chains)
    angles = sizehint!(Float64[0, 180], chains .|> length |> sum)

    for chain in chains
        chainConVectors = chain[2:end] .- chain[1:end-1] .|> normalize
        chainAngles = dot.(-chainConVectors[2:end], chainConVectors[1:end-1]) .|> acos .|> rad2deg

        append!(angles, chainAngles)
    end

    println()
    histogram(angles, vertical = true, title="Bond Angle Distribution", nbins=90, border= :bold) |> show
    println("\n")
end