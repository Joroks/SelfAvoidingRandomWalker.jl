"""
    randomWalk(chainLengths, boxSize, bondLength, minDistance, targetAngle, maxNumTries, particles = [])

Run the self-avoiding random walk algorithm in a periodic simulation box.

# Arguments
- `chainLenghts`: the lenghts of the chains that will be placed within the simulation box.
- `boxSize`: the size of the simulation box specified by a 3x3 Matrix. The columns of the matrix specify the basis vectors spanning the parallelepiped representing the simulation box.
- `bondLenght`: the distance between to consecutive beads of a chain
- `minDistance`: smallest permitted distance between any two beads. The previous bead in a chain is ignored.
- `targetAnlge`: target for the resulting angle enclosed with the two previous beads. This can also be a range `(targetMin, targetMax)` and can also include a range of allowed angles `(targetMin, targetMax, allowedMin, allowedMax)`
- `particles`: a list of particles that will be avoided by the random walker
- `maxNumTries`: specifies how often the algorithm tries to place any given bead before backtracking one step.
"""
function randomWalk(chainLengths, boxSize, bondLength, minDistance, targetAngle=(0, 180), particles = Particle[], maxNumTries=50)
    boxSize = Matrix3x3(boxSize)

    numAtoms = sum(chainLengths)

    sort!(chainLengths, by=-)

    chains = sizehint!(Vector{Vector3}[], length(chainLengths))

    cellGrid = CellGrid(boxSize, minDistance)

    initialRD = searchPattern(0, 180)
    growthRD = searchPattern(targetAngle...)

    desc = "Running random walker:"
    barlen = ProgressMeter.tty_width(desc, stderr, true) - 1 # workaround for incorrect length calculation
    p = Progress(numAtoms, desc, showspeed=true, barlen=barlen)

    numPlaced = 0
    highestNumPlaced = 0
    numTries = 0

    for targetLength in chainLengths
        points = sizehint!(Vector3[], targetLength)
        directions = sizehint!(Matrix3x3[], targetLength)
        tries = sizehint!(Int64[], targetLength)

        while length(points) < targetLength

            while !isempty(points) && tries[end] >= maxNumTries
                remove!(cellGrid, points[end])

                pop!(points)
                pop!(directions)
                pop!(tries)
                numPlaced -= 1
            end

            if isempty(points)
                nextStep = boxSize*rand(3)
                nextDirection = initialRD(0)
            else
                tryFactor = (tries[end] + rand()) / maxNumTries
                nextDirection = directions[end] * growthRD(tryFactor)
                nextStep = points[end] + nextDirection[:,1] * bondLength

                tries[end] += 1
            end

            ignore(n) = !isempty(points) && n == points[end]

            if !hasNeighbours(cellGrid, nextStep, ignore) && !any(particles .|> p -> isInParticle(nextStep, boxSize, p))
                add!(cellGrid, nextStep)

                push!(points, nextStep)
                push!(directions, nextDirection)
                push!(tries, 0)
                highestNumPlaced = max(numPlaced += 1, highestNumPlaced)
            end

            numTries += 1
            update!(p, highestNumPlaced, showvalues = [("atoms placed", highestNumPlaced), ("number of tries", numTries), ("average number of tries", numTries/highestNumPlaced)])
        end

        push!(chains, points)
    end

    finish!(p)
    return chains
end