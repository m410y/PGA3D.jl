export AbstractMVec, AbstractEVec

abstract type AbstractMVec{T<:Real} <: Number end
abstract type AbstractEVec{T<:Real} <: AbstractMVec{T} end

@assert length(indices) == D
@assert length(basis) == exp2(D)

function grade_basis(grade)
    basis[findall(sym -> length(String(sym)) - 1 in grade, basis)]
end

# constructors
for (type, grade) in type_grade
    @eval export $type
    @eval struct $type{T<:Real} <: $(all(iseven.(grade)) ? AbstractEVec : AbstractMVec){T}
        $(map(sym -> :($sym::T), grade_basis(grade))...)
    end
    @eval $type(m::$type) = m
    @eval $type{T}(m::$type) where {T<:Real} = $type{T}(m)
end

# promotion
for (type, grade) in type_grade
    if 0 in grade
        @eval $type(x::T) where {T<:Real} = $type{T}(x)
        @eval @generated $type{T}(x::Real) where {T<:Real} = Expr(:call, $type{T},
            map(fieldnames($type)) do field
                field == Symbol(basis_char) ? :(x) : 0
        end...)
    end
    for (type2, grade2) in type_grade
        if grade2 in grade && type2 != type
            @eval $type(m::$type2{T}) where {T} = $type{T}(m::$type2)
            @eval @generated $type{T}(m::$type2) where {T<:Real} = Expr(:call, $type{T},
                map(fieldnames($type)) do field
                    hasfield($type2, field) ? :(m.$field) : 0
            end...)
        end
    end
end
function Base.promote_rule(
    ::Type{E},
    ::Type{S},
) where {E<:AbstractEVec{T},S<:Real} where {T<:Real}
    EVec{promote_type(T, S)}
end
function Base.promote_rule(
    ::Type{E1},
    ::Type{E2},
) where {E1<:AbstractEVec{T},E2<:AbstractEVec{S}} where {T<:Real,S<:Real}
    EVec{promote_type(T, S)}
end
function Base.promote_rule(
    ::Type{M},
    ::Type{S},
) where {M<:AbstractMVec{T},S<:Real} where {T<:Real}
    MVec{promote_type(T, S)}
end
function Base.promote_rule(
    ::Type{M1},
    ::Type{M2},
) where {M1<:AbstractMVec{T},M2<:AbstractMVec{S}} where {T<:Real,S<:Real}
    MVec{promote_type(T, S)}
end

# scalar projection
Base.real(m::T) where {T<:AbstractMVec} = hasfield(T, Symbol(basis_char)) ? m.e : 0
Base.real(M::Type{T}) where {T<:AbstractMVec} =
    hasfield(M, Symbol(basis_char)) ? fieldtype(M, Symbol(basis_char)) :
    error("$T has no real part")

# boolean
Base.isreal(m::T) where {T<:AbstractMVec} =
    all((name != Symbol(basis_char) && iszero(getfield(m, name)) for name in fieldnames(T)))
Base.isinteger(m::AbstractMVec{T}) where {T} = isreal(m) & T <: Integer
Base.isfinite(m::T) where {T<:AbstractMVec} =
    all(isfinite(getfield(m, name)) for name in fieldnames(T))
Base.isnan(m::T) where {T<:AbstractMVec} =
    any(isnan(getfield(m, name)) for name in fieldnames(T))
Base.isinf(m::T) where {T<:AbstractMVec} =
    all(isinf(getfield(m, name)) for name in fieldnames(T))
Base.iszero(m::T) where {T<:AbstractMVec} =
    all(iszero(getfield(m, name)) for name in fieldnames(T))
Base.isone(m::T) where {T<:AbstractMVec} = all((
    name == Symbol(basis_char) ? isone(m.e) : iszero(getfield(m, name)) for
    name in fieldnames(T)
))
==(m1::T, m2::T) where {T<:AbstractMVec} =
    all((getfield(m1, name) == getfield(m2, name) for name in fieldnames(T)))
Base.in(m::AbstractMVec, r::AbstractRange{<:Real}) = isreal(m) && real(m) in r

# utility
Base.flipsign(m::AbstractMVec, x::Real) = ifelse(signbit(x), -m, m)
Base.bswap(m::T) where {T<:AbstractMVec} =
    T((bswap(getfield(m, name)) for name in fieldnames(T))...)
for type in keys(type_grade)
    @eval Base.widen(::Type{$type{T}}) where {T} = $type{widen(T)}
    @eval Base.float(::Type{$type{T}}) where {T<:AbstractFloat} = $type{T}
    @eval Base.float(::Type{$type{T}}) where {T} = $type{float(T)}
end
