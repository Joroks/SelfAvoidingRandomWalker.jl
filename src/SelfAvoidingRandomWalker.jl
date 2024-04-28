module SelfAvoidingRandomWalker

using LinearAlgebra
using StaticArrays
using ProgressMeter
using Distributions
using UnicodePlots
using Printf

const Vector3 = SVector{3, Float64}
const Matrix3x3 = SMatrix{3, 3, Float64}

include("RandomWalker.jl")
include("CellGrid.jl")
include("SearchPattern.jl")
include("ChainLengths.jl")
include("SaveAsLammps.jl")
include("Particle.jl")

export randomWalk
export saveAsLammps
export randomLengths
export Particle
export printChainLengthStatistics
export printAngleStatistics

end # module SARW