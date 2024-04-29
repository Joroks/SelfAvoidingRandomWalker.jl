# SelfAvoidingRandomWalker

[![Build Status](https://github.com/Jorojs/SelfAvoidingRandomWalker.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Jorojs/SelfAvoidingRandomWalker.jl/actions/workflows/CI.yml?query=branch%3Amaster)

SelfAvoidingRandomWalker is a julia package for quickly generating initial datasets for use in Molecular Dynamics simulations of linear Thermoplastic Polymers in a periodic simulation box. 


## Installation

Within julia, execute

```julia
using Pkg; Pkg.add(url="https://github.com/Joroks/SelfAvoidingRandomWalker.jl")
```



## Model parameters

- `boxSize`: the size of the simulation box specified by a 3x3 Matrix. The columns of the matrix specify the basis vectors spanning the parallelepiped representing the simulation box.
It has to have the following shape in order to be compatible with the Lammps data format:

```
    xhi  xy   xz
    0    yhi  yz
    0    0    zhi
```
with
```
    xhi, yhi, zhi > 0

    -xhi/2 ≤ xy ≤ xhi/2
    -xhi/2 ≤ xz ≤ yhi/2
    -xhi/2 ≤ yz ≤ yhi/2
```

- `bondLenght`: the distance between to consecutive beads of a chain
- `minDistance`: smallest permitted distance between any two beads. The previous bead in a chain is ignored.
- `targetAnlge`: target for the resulting angle enclosed with the two previous beads. This can also be a range `(targetMin, targetMax)` and can also include a range of allowed angles `(targetMin, targetMax, allowedMin, allowedMax)`
- `numAtoms`: the number of Atoms to be placed within the simulation box
- `averageLength`: the target for the average length of the chains to be placed within the simulation box
- `dispersity`: the target for the dispersity of the chains to be placed within the simulation box
- `particles`: a list of particles that will be avoided by the random walker. A particle is described by three fields:
    - `position` (3D-Vector): the position of the center of the particle
    - `shape` (3x3-Matrix): in most cases this will be a diagonal matrix where each entry describes the radius of that particle in the corresponding direction
    - `pvalue` (3D-Vector): the exponent used for the calculation of distance between a point and the particle
        - the following equation is used to determine wheter a point is valid : `isvalid(p) = sum(shape\(p-position) .^ pvalue) > 1`
        - in reality the impementation is slightly different to account for periodicity
- `maxNumTries`: specifies how often the algorithm tries to place any given bead before backtracking one step.


## Usage

```julia
using SelfAvoidingRandomWalker
using LinearAlgebra # for diagm()

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

## Generating Chain Lenghts

chainLengths = randomLengths(numAtoms, averageLength, dispersity)
printChainLengthStatistics(chainLengths)

## Running the SARW algorithm

points = randomWalk(
    chainLengths,
    boxSize,
    bondLength,
    minDistance,
    targetAngle,
    particles,
    maxNumTries
)


printAngleStatistics(points)

## Save result

saveAsLammps(boxSize, points, "<output-filename>", "#<Lammps-Header>")
```