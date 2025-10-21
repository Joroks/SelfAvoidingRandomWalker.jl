module SelfAvoidingRandomWalker

using LinearAlgebra
using StaticArrays
using ProgressMeter
using Distributions
using UnicodePlots
using Printf
using Pipe: @pipe
import ModernGL, GLFW
import CImGui as ig

const Vector3 = SVector{3, Float64}
const Matrix3x3 = SMatrix{3, 3, Float64, 9}

function visualize_search_pattern end

include("CountingBins.jl")
include("RandomWalker.jl")
include("SearchPattern.jl")
include("ChainLengths.jl")
include("SaveAsLammps.jl")

export SARW, run_SARW!, Particle, randomLengths,  printChainLengthStatistics, printAngleStatistics, saveAsLammps, chain_generator

end # module SelfAvoidingRandomWalker