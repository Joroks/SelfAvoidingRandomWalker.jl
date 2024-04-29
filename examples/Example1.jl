using LinearAlgebra
using SelfAvoidingRandomWalker

## model Parameters
boxSize = [
    200 10 20
    0 130 30
    0 0 100
]

minDistance = 0.9
bondLength = 0.7
targetAngle = (160, 160)

numAtoms = 1_120_042
averageLength = 1000
dispersity = 2.4

particles = Particle[
    Particle([30, 0, 0], 15) # sphere
    Particle([100, 50, 10], diagm([30, 10, 10])) # ellipse
    Particle(boxSize*[0.5, 0.5, 0.5], 10, Inf) # cube
    Particle([0, 30, 0], diagm([30, 4, 30]), [2, Inf, 2]) # flat round disc
    Particle(boxSize*[0.25, 0.25, 0.25], diagm([5, 5, 40]), [2, 2, Inf]) # cylinder
]

maxNumTries = 50

## executing SARW algorithm - do not edit
chainLengths = randomLengths(numAtoms, averageLength, dispersity)
printChainLengthStatistics(chainLengths)

points = randomWalk(chainLengths, boxSize, bondLength, minDistance, targetAngle, particles, maxNumTries)
printAngleStatistics(points)

if !isempty(ARGS)
    saveAsLammps(boxSize, points, "output/test.data", "#")
end