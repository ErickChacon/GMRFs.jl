"""
    AbstractGMRF

Sypertype of Gaussian Markov random fields with mean `μ` and precison matrix `Q = κS`,
where `S` is the structure matrix and `κ` the scale parameter.
"""
abstract type AbstractGMRF <: Distributions.ContinuousMultivariateDistribution end

"""
    length(d::AbstractGMRF)

Return the sampling dimension of the GMRF `d`.
"""
Base.length(d::AbstractGMRF)

"""
    scale(d::AbstractGMRF)

Return the scale parameter of the GMRF `d`.
"""
Distributions.scale(d::AbstractGMRF)

"""
    structure(d::AbstractGMRF)

Returns the structure matrix of the GMRF `d`.
"""
function structure end

precision(d::AbstractGMRF) = scale(d) * structure(d)


## Random generator

function Distributions._rand!(rng::AbstractRNG, d::AbstractGMRF, x::AbstractVector{T}) where T<:Real
    randn!(rng, x)
    copyto!(x, cholesky(structure(d)).UP \ x)
    ldiv!(sqrt(scale(d)), x)
    return x
end

function Distributions._rand!(rng::AbstractRNG, d::AbstractGMRF, x::AbstractArray{<:Real})
    U = cholesky(structure(d)).UP
    @inbounds for xi in Distributions.eachvariate(x, Distributions.variate_form(typeof(d)))
        randn!(rng, xi)
        copyto!(xi, U \ xi)
        ldiv!(sqrt(scale(d)), xi)
    end
    return x
end

# This functions is required for the default Distributions.rand!. This is a left-division
# method for a sparse cholesky facor and a subarray.
function LinearAlgebra.:\(L::SuiteSparse.CHOLMOD.FactorComponent,
                          b::SuiteSparse.CHOLMOD.AbstractVector)
    reshape(Matrix(L \ SuiteSparse.CHOLMOD.Dense(b)), length(b))
end

## Logarithm of the pdf

function Distributions._logpdf(d::AbstractGMRF, x::AbstractVector{<:Real})
    n = length(d)
    chol = cholesky(structure(d))
    logpdf = -0.5 * n * log(2.0 * pi)
    logpdf += 0.5 * (n * log(scale(d)) + LinearAlgebra.logdet(chol))
    logpdf -= 0.5 * scale(d) * x' * structure(d) * x
    return(logpdf)
end

# TODO: Find a way to not repeat this function from Distributions
@inline function Distributions.logpdf(d::AbstractGMRF, x::AbstractArray{<:Real,1})
    @boundscheck begin
        size(x) == size(d) ||
            throw(DimensionMismatch("inconsistent array dimensions"))
    end
    return Distributions._logpdf(d, x)
end

# TODO: Custom modification from to Distributions.logpdf
@inline function Distributions.logpdf(d::AbstractGMRF, x::AbstractArray{<:Real,M}) where {M}
    @boundscheck begin
        M > 1 ||
            throw(DimensionMismatch(
                "number of dimensions of `x` ($M) must be greater than number of dimensions of `d` (1)"
            ))
        ntuple(i -> size(x, i), Val(1)) == size(d) ||
            throw(DimensionMismatch("inconsistent array dimensions"))
    end

    n = length(d)
    chol = cholesky(structure(d))
    logpdf = -0.5 * n * log(2.0 * pi)
    logpdf += 0.5 * (n * log(scale(d)) + LinearAlgebra.logdet(chol))
    lpdf = @inbounds map(xi -> -0.5 * scale(d) * xi' * structure(d) * xi,
               Distributions.eachvariate(x, Distributions.variate_form(typeof(d))))
    return logpdf .+ lpdf
end

## IO

Base.show(io::IO, d::AbstractGMRF) = print(io, "$(length(d)) GMRF")

function Base.show(io::IO, ::MIME"text/plain", d::AbstractGMRF)
    println(io, d)
    println(io, "  S: ", summary(structure(d)))
    print(io, "  κ: ", scale(d))
end

## GMRF implementations

for filename in ["gmrf.jl", "cgmrf.jl"]
    include(joinpath("gmrf", filename))
end

