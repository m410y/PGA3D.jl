export plane, point, direction

plane(a::Real, b::Real, c::Real, d::Real) = Vec(a, b, c, d)
point(x::Real, y::Real, z::Real, w::Real = 1) = PVec(x, y, z, w)
function point(p::AbstractVector)
    @assert length(p) == 3
    point(p..., 1)
end
direction(x::Real, y::Real, z::Real) = point(x, y, z, 0)
function direction(v::AbstractVector)
    @assert length(v) == 3
    point(v..., 0)
end
line(vx::Real, vy::Real, vz::Real, px::Real, py::Real, pz::Real) = BiVec(vx, vy, vz, px, py, pz)
function line(v::AbstractVector, p::AbstractVector)
    @assert length(v) == length(p) == 3
    BiVec(v..., p...)
end
axis(vx::Real, vy::Real, vz::Real) = line(vx, vy, vz, 0, 0, 0)
