using PolygonOps
using GeophysicalModelGenerator
using Base.Threads


function AddPoly!(Phase, Temp, Grid::AbstractGeneralGrid;
    vertices::AbstractVector,
    phase=ConstantPhase(1),
    T=nothing)
    @assert length(vertices) >= 3 "At least 3 vertices are required to form a polygon."

    # Determine dimensionality of vertices
    dim = length(vertices[1])
    @assert dim == 2 || dim == 3 "Vertices must have 2 or 3 coordinates."

    # Retrieve 3D data arrays for the grid
    X, Y, Z = coordinate_grids(Grid)



    polygon = [[vertex[1], vertex[2]] for vertex in vertices]
    points = [(x, z) for (x, z) in zip(X, Z)]
    inside = Vector{Bool}(undef, length(points))
    @inbounds @threads for i in eachindex(points)
        inside[i] = inpolygon(points[i], polygon, in = 1, on = 1, out = 0)
    end
    ind = findall(inside)

    if isempty(ind)
        # println("points: ", points)
        println("No points are inside the polygon with vertices: ", vertices)
        return nothing
    end
    # Compute thermal structure accordingly
    if T !== nothing && !isempty(ind)
        Temp[ind] = Compute_ThermalStructure(Temp[ind], X[ind], Y[ind], Z[ind], T)
    end

    # Set the phase
    if !isempty(ind)
        Phase[ind] = Compute_Phase(Phase[ind], Temp[ind], X[ind], Y[ind], Z[ind], phase)
    end

    return nothing
end

AddPoly!(model::Model; kwargs...) = AddPoly!(model.Grid.Phases, model.Grid.Temp, model.Grid.Grid; kwargs...)
