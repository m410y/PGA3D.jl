module PGA3D

import Base: +, -, *, ~, /, ==
using LinearAlgebra
import LinearAlgebra: ⋅, ×

const metric = (0, 1, 1, 1)
const D = length(metric)
const basis_char = 'e'
const indices = '0':'3'
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
const type_grade = Dict(
    :Vec => 1,
    :BiVec => 2,
    :PVec => D - 1,
    :PScalar => D,
    :EVec => 0:2:D,
    :MVec => 0:D,
)

include("core.jl")
include("generic.jl")
include("handwrite.jl")
include("primitives.jl")

end # module PGA3D
