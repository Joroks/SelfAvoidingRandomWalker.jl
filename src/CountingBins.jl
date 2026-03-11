struct Bin{P}
    parent::P
    index::Int
end

Base.IteratorSize(::Bin) = Base.SizeUnknown()
Base.eltype(::Bin) = Int
function Base.iterate(self::Bin, state=nothing)
    if state === nothing
        @inbounds state = self.parent.bins[self.index]
    end
    
    state < 0 && return nothing
    state, @inbounds self.parent.entries[state][1]
end

function Base.show(io::IO, self::Bin)
    Base.show_vector(io, collect(self))
end

mutable struct CountingBins{N} <: AbstractArray{Bin, N}
    const bins::Memory{Int}
    const entries::Memory{Int}
    const size::NTuple{N, Int}
    counter::Int

    function CountingBins(capacity::Int, dims::NTuple{N, Int}) where N
        self = new{N}(
            Memory{Int}(undef, prod(dims)),
            Memory{Int}(undef, capacity),
            dims, 0
        )

        self.bins .= -(1:prod(dims))
        return self
    end
end

Base.size(self::CountingBins) = self.size
Base.IndexStyle(::CountingBins) = Base.IndexLinear()

function Base.getindex(self::CountingBins, i::Int)
    @boundscheck checkbounds(self, i)
    return Bin(self, i)
end

function increment!(self::Bin)
    count = self.parent.counter += 1
    self.parent.entries[count] = self.parent.bins[self.index]
    self.parent.bins[self.index] = count
    return count
end

function decrement!(self::CountingBins)
    next = self.entries[self.counter]
    bin = next
    while bin > 0; bin = self.entries[bin] end
    self.bins[-bin] = next
    self.counter -= 1
end

Base.count(self::CountingBins) = self.counter