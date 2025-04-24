export plane, point, direction, rotaxis

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
rotaxis(vx::Real, vy::Real, vz::Real, px::Real, py::Real, pz::Real) = normalize(direction(vx, vy, vz)) ∨ point(px, py, pz)
rotaxis(vx::Real, vy::Real, vz::Real) = normalize(direction(vx, vy, vz)) ∨ point(0, 0, 0)
rotaxis(v::AbstractVector, p::AbstractVector) = normalize(direction(v)) ∨ point(p)
rotaxis(v::AbstractVector) = normalize(direction(v)) ∨ point(0, 0, 0)
