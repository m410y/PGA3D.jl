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

for type in keys(type_grade)
    @eval @generated *(m::$type, x::Real) = Expr(:call, $type, map(fieldnames($type)) do field
        :(m.$field * x)
    end...)
    @eval @generated *(x::Real, m::$type) = Expr(:call, $type, map(fieldnames($type)) do field
        :(x * m.$field)
    end...)
    @eval @generated /(m::$type, x::Real) = Expr(:call, $type, map(fieldnames($type)) do field
        :(m.$field / x)
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
