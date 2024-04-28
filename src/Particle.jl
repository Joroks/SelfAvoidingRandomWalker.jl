struct Particle
    position::SVector{3, Float64}
    shape::SMatrix{3,3, Float64}
    pvalue::SVector{3, Float64}

    Particle(position::AbstractVector, radius::Real) = new(position, radius*I(3), [2, 2, 2])
    Particle(position::AbstractVector, shape::AbstractArray) = new(position, shape, [2, 2, 2])
    Particle(position::AbstractVector, shape::AbstractArray, pvalue::Real) = new(position, shape, ones(3)*pvalue)
    Particle(position::AbstractVector, radius::Number, pvalue::Real) = new(position, radius*I(3), ones(3)*pvalue)
    Particle(position::AbstractVector, shape::AbstractArray, pvalue::AbstractVector) = new(position, shape, pvalue)
end

isInParticle(point::AbstractVector, boxSize::AbstractMatrix, particle::Particle) = particle.shape\roundMod(point - particle.position, boxSize) |> x -> abs.(x).^particle.pvalue |> sum |> x -> x <= 1
