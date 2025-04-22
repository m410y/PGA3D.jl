export AbstractMVec, AbstractEVec
export ⋆, ⋅, ∧, ×, ∨

abstract type AbstractMVec{T<:Real} <: Number end
abstract type AbstractEVec{T<:Real} <: AbstractMVec{T} end

const metric = (0, 1, 1, 1)
const D = length(metric)
const basis_char = 'e'
const indices = '0':'3'
@assert length(indices) == D
const basis = [
    :e,
    :e1,
    :e2,
    :e3,
    :e0,
    :e23,
    :e31,
    :e12,
    :e01,
    :e02,
    :e03,
    :e032,
    :e013,
    :e021,
    :e123,
    :e0123,
]
@assert length(basis) == exp2(D)
const type_grade = Dict(
    :Vec => 1,
    :BiVec => 2,
    :PVec => D - 1,
    :PScalar => D,
    :EVec => 0:2:D,
    :MVec => 0:D,
)

function grade_basis(grade)
    basis[findall(sym -> length(String(sym)) - 1 in grade, basis)]
end

for (type, grade) in type_grade
    @eval export $type
    @eval struct $type{T<:Real} <: $(all(iseven.(grade)) ? AbstractEVec : AbstractMVec){T}
        $(map(sym -> :($sym::T), grade_basis(grade))...)
    end

    if length(grade_basis(grade)) > 1
        @eval $type(args::Vararg{Real,fieldcount($type)}) =
            $type{promote_type(typeof.(args)...)}(args...)
        if 0 in grade
            @eval $type(x::T) where {T<:Real} = $type{T}(x)
            @eval $type{T}(x::Real) where {T<:Real} = $type{T}(
                (
                    field == Symbol(basis_char) ? x : zero(T) for
                    field in fieldnames($type)
                )...,
            )
        end
    end
    @eval $type(m::AbstractMVec{T}) where {T} = $type{T}(m)
    @eval $type{T}(m::M) where {T<:Real,M<:AbstractMVec} = $type{T}(
        (
            hasfield(M, field) ? getfield(m, field) : zero(T) for
            field in fieldnames($type)
        )...,
    )

    @eval Base.widen(::Type{$type{T}}) where {T} = $type{widen(T)}
    @eval Base.float(::Type{$type{T}}) where {T<:AbstractFloat} = $type{T}
    @eval Base.float(::Type{$type{T}}) where {T} = $type{float(T)}
end

# promotion
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

# addition and subtraction
for type in keys(type_grade), op in (:+, :-)
    @eval @generated $op(m::$type) = Expr(:call, $type, map(fieldnames($type)) do field
        :($$op(m.$field))
    end...)
    @eval @generated $op(m1::$type, m2::$type) =
        Expr(:call, $type, map(fieldnames($type)) do field
            :($$op(m1.$field, m2.$field))
        end...)
end

function _basisreduce(arr::Vector{Char})
    count = zeros(Int8, D)
    if isempty(arr)
        return 1, Symbol(basis_char), Tuple(count)
    end
    # insertion sort
    count[1+first(arr)-first(indices)] += 1
    perms = 0
    for i = 2:length(arr)
        j = i
        count[1+arr[j]-first(indices)] += 1
        while j > 1 && arr[j-1] > arr[j]
            arr[j-1], arr[j] = arr[j], arr[j-1]
            perms += 1
            j -= 1
        end
    end
    res = [basis_char]
    sizehint!(res, D + 1)
    for (n, c) in zip(count, indices)
        if isodd(n)
            push!(res, c)
        end
    end
    sign = iseven(perms) ? 1 : -1
    sign, Symbol(String(res)), Tuple(count .÷ 2)
end
_basisreduce(a::Symbol) = _basisreduce(collect(String(a)[2:end]))
_basisreduce(a::Symbol, b::Symbol) = _basisreduce(
    collect(String(a)[2:end] * String(b)[2:end]),
)
_grade(a::Symbol) = length(unique(String(a))) - 1
_dual(a::Symbol) = _basisreduce(a, only(grade_basis(D)))[1:2]

function _grade_type(grade)
    types = findall(type_grade) do _grade
        isempty(setdiff(grade, _grade))
    end
    argmin(types) do type
        length(grade_basis(type_grade[type]))
    end
end

function _unary_expr(basisop, type::Type{<:AbstractMVec})
    optable = Dict()
    grade_set = Set()
    for name in fieldnames(type)
        sign, sym = basisop(name)
        if sign != 0
            optable[sym] = sign, name
            push!(grade_set, _grade(sym))
        end
    end
    if length(grade_set) == 0
        return 0
    end
    if length(grade_set) == 1 && only(grade_set) == 0
        sign, sym = optable[Symbol(basis_char)]
        return sign > 0 ? :(m.$sym) : :(-m.$sym)
    end
    res_type = _grade_type(sort(collect(grade_set)))
    Expr(:call, res_type, map(grade_basis(type_grade[res_type])) do _name
        sign1, name, _ = _basisreduce(_name)
        sign2, sym = optable[name]
        sign = sign1 * sign2
        sign > 0 ? :(m.$sym) : :(-m.$sym)
    end...)
end

# hodge star
⋆(x::Real) = PScalar(x)
for type in keys(type_grade)
    @eval @generated function ⋆(m::$type)
        _unary_expr(_dual, $type)
    end
end

# reverse
for type in keys(type_grade)
    @eval @generated function ~(m::$type)
        Expr(:call, $type, map(fieldnames($type)) do sym
            g = _grade(sym)
            iseven((g * (g - 1)) ÷ 2) ? :(m.$sym) : :(-m.$sym)
        end...)
    end
end

function _mul_terms(pos, neg)
    if isempty(pos) && isempty(neg)
        return 0
    end
    posex, negex = map((pos, neg)) do terms
        if isempty(terms)
            return nothing
        end
        termex = map(terms) do term
            :(m1.$(term[1]) * m2.$(term[2]))
        end
        if length(termex) == 1
            only(termex)
        else
            Expr(:call, :+, termex...)
        end
    end
    if isnothing(negex)
        :($posex)
    elseif isnothing(posex)
        :(-$negex)
    else
        :($posex - $negex)
    end
end

function _mul_expr(basisprod, t1::Type{<:AbstractMVec}, t2::Type{<:AbstractMVec})
    multable = Base.ImmutableDict(
        [_basisreduce(name)[2] => (pos = [], neg = []) for name in basis]...,
    )
    grade_set = Set()
    for n1 in fieldnames(t1), n2 in fieldnames(t2)
        sign, sym = basisprod(n1, n2)
        if sign > 0
            push!(multable[sym].pos, (n1, n2))
            push!(grade_set, _grade(sym))
        elseif sign < 0
            push!(multable[sym].neg, (n1, n2))
            push!(grade_set, _grade(sym))
        end
    end
    if length(grade_set) == 0
        return 0
    end
    if length(grade_set) == 1 && only(grade_set) == 0
        return _mul_terms(multable[Symbol(basis_char)]...)
    end
    type = _grade_type(sort(collect(grade_set)))
    Expr(:call, type, map(grade_basis(type_grade[type])) do _name
        sign, name, _ = _basisreduce(_name)
        if sign > 0
            _mul_terms(multable[name]...)
        else
            _mul_terms(reverse(multable[name])...)
        end
    end...)
end

function _geomprod(a::Symbol, b::Symbol)
    sign, res, count = _basisreduce(a, b)
    sign *= prod(metric .^ count)
    sign, res
end

for type1 in keys(type_grade), type2 in keys(type_grade)
    @eval @generated function *(m1::$type1, m2::$type2)
        _mul_expr(_geomprod, $type1, $type2)
    end
    @eval @generated function ⋅(m1::$type1, m2::$type2)
        _mul_expr($type1, $type2) do a, b
            sign, res = _geomprod(a, b)
            sign *= _grade(res) == abs(_grade(a) - _grade(b))
            sign, res
        end
    end
    @eval @generated function ∧(m1::$type1, m2::$type2)
        _mul_expr($type1, $type2) do a, b
            sign, res = _geomprod(a, b)
            sign *= _grade(res) == _grade(a) + _grade(b)
            sign, res
        end
    end
    @eval @generated function ×(m1::$type1, m2::$type2)
        _mul_expr($type1, $type2) do a, b
            sign1, res = _geomprod(a, b)
            sign2, _ = _geomprod(b, a)
            sign = (sign1 - sign2) ÷ 2
            sign, res
        end
    end
    @eval @generated function ∨(m1::$type1, m2::$type2)
        _mul_expr($type1, $type2) do a, b
            sign1, da = _dual(a)
            sign2, db = _dual(b)
            sign3, dres = _geomprod(da, db)
            sign3 *= _grade(dres) == _grade(da) + _grade(db)
            sign4, res = _dual(dres)
            sign1 * sign2 * sign3 * sign4, res
        end
    end
end

# norm and inverse
Base.abs2(m::AbstractMVec) = real(m ⋅ ~m)
Base.abs(m::AbstractMVec) = sqrt(abs(abs2(m)))
LinearAlgebra.norm(m::AbstractMVec) = abs(m)
LinearAlgebra.normalize(m::AbstractMVec) = iszero(norm(m)) ? m : m / norm(m)
Base.inv(m::AbstractMVec) = ~m / abs2(m)
/(x::Real, m::AbstractMVec) = x * inv(m)
/(m1::AbstractMVec, m2::AbstractMVec) = m1 * ~m2 / abs2(m2)

# exponentiation
function Base.exp(v::Vec{T}) where {T}
    s = sqrt(abs2(v))
    iszero(s) ? one(T) : cosh(s) + v * sinh(s) / s
end
function Base.exp(B::BiVec{T}) where {T}
    s = sqrt(-abs2(B))
    p = B ∧ B / 2
    cos(s) + B * sinc(s / π) - B * p * cosc(s / π) / π + p * sinc(s / π)
end
function Base.exp(V::PVec{T}) where {T}
    s = abs(V.e123)
    iszero(s) ? one(T) : cos(s) + V * sinc(s / π)
end
Base.exp(S::PScalar{T}) where {T} = one(T) + S
Base.exp(m::EVec{T}) where {T} = exp(real(m)) * exp(BiVec{T}(m)) * exp(PScalar{T}(m))
