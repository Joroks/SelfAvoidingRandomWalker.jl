struct CellGrid
    boxSize::Matrix3x3
    cellSize::Matrix3x3
    numCells::Matrix3x3

    searchRadius::Float64
    overlapMask::Vector{Vector3}

    cells::Dict{Vector3, Vector{Vector3}}
end

function CellGrid(boxSize, searchRadius)
    numCells = round.(boxSize/searchRadius, RoundFromZero)
    cellSize = boxSize/numCells
    overlapMask = gridVectors([-1 0 1], 3) |> eachcol .|> Vector3
    cells = Dict{Vector3, Vector{Vector3}}()

    return CellGrid(boxSize, cellSize, numCells, searchRadius, overlapMask, cells)
end

const emptyVecArray = Vector3[]

function hasNeighbours(grid::CellGrid, point, ignore = x -> false)
    index = cellIndex(point, grid.cellSize)

    isInRadius(p) = periodicDistance(p, point, grid.boxSize) <= grid.searchRadius    

    Iterators.map(o -> floorMod(index + o, grid.numCells), grid.overlapMask) |>
        indices -> Iterators.flatmap(i -> get(grid.cells, i, emptyVecArray), indices) |>
        points -> Iterators.map(p -> isInRadius(p) && !ignore(p) , points) |>
        any
end

function add!(grid::CellGrid, point::Vector3)
    index = cellIndex(point, grid.cellSize)
    indexMod = floorMod(index, grid.numCells)

    cell = get!(grid.cells, indexMod, Vector3[])

    push!(cell, point)

    return grid
end

function remove!(grid::CellGrid, point)
    index = floorMod(cellIndex(point, grid.cellSize), grid.numCells)

    cell = get(grid.cells, index, emptyVecArray)

    firstOccurance = findfirst(p -> p==point, cell)
    deleteat!(cell, firstOccurance)

    return grid
end

function gridVectors(a, n)
    if n == 1
        return a'
    end

    npart = length(a)
    partlength = npart^(n-1)
    ncols = partlength*npart

    prevBlock = gridVectors(a, n-1)
    currBlock = similar(a, n, ncols)

    for part in 1:npart
        partIndices = (1:partlength) .+ partlength*(part-1)

        currBlock[1:(n-1), partIndices] = prevBlock
        currBlock[n, partIndices] .= a[part]
    end

    return currBlock
end

cellIndex(point, cellSize) = floor.(cellSize \ point)

floorMod(index, numCells) = index - numCells * floor.(numCells \ index)

roundMod(point, boxSize) = point - boxSize * round.(boxSize \ point)

periodicDistance(a, b, boxSize) = roundMod(b - a, boxSize) |> norm
