export >>>

# norm and inverse
Base.abs2(m::AbstractMVec) = real(m ⋅ ~m)
Base.abs(m::AbstractMVec) = sqrt(abs(abs2(m)))
LinearAlgebra.norm(m::AbstractMVec) = abs(m)
LinearAlgebra.normalize(m::AbstractMVec) = iszero(abs(m)) ? m : m / abs(m)
Base.inv(m::AbstractMVec) = ~m / abs2(m)
/(x::Real, m::AbstractMVec) = x * inv(m)
/(m1::AbstractMVec, m2::AbstractMVec) = m1 * ~m2 / abs2(m2)

# exponentiation
function Base.exp(v::Vec{T}) where {T}
    s = abs(v)
    iszero(s) ? one(T) : cosh(s) + v * sinh(s) / s
end
function Base.exp(B::BiVec{T}) where {T}
    s = abs(B)
    p = B ∧ B / 2
    cos(s) + B * sinc(s / π) - B * p * cosc(s / π) / π + p * sinc(s / π)
end
function Base.exp(V::PVec{T}) where {T}
    s = abs(V.e123)
    iszero(s) ? one(T) : cos(s) + V * sinc(s / π)
end
Base.exp(S::PScalar{T}) where {T} = one(T) + S
Base.exp(m::EVec{T}) where {T} = exp(real(m)) * exp(BiVec{T}(m)) * exp(PScalar{T}(m))

>>>(m::AbstractEVec, v::T) where {T<:AbstractMVec} = T(m * v * ~m)
