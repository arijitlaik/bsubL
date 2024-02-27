using Parameters

# Define the abstract type for 2D shapes
abstract type Shape2D end

# Define the unique fields for Slab2D
@with_kw mutable struct Slab2D <: Shape2D
    top_x::Float64
    top_y::Float64
    length::Float64
    taper::Float64
    dip::Float64
    depth::Float64
    thickness_array::Vector{Float64}
    polygons::Union{Nothing, Vector{Vector{Vector{Float64}}}} = nothing
end

# Define the unique fields for Indentor2D
@with_kw mutable struct Indentor2D <: Shape2D
    top_x::Float64
    top_y::Float64
    length::Float64
    taper::Float64
    taper2::Union{Nothing, Float64} = nothing
    thickness_array::Vector{Float64}
    polygons::Union{Nothing, Vector{Vector{Vector{Float64}}}} = nothing
end

@with_kw mutable struct OverRidingPlate2D <: Shape2D
    top_x::Float64
    top_y::Float64
    length::Float64
    taper::Float64
    dip::Float64
    thickness_array::Vector{Float64}
    polygons::Union{Nothing, Vector{Vector{Vector{Float64}}}} = nothing
end

# Define a function to compute vertices for Slab2D
function compute_poly(slab::Slab2D)
    @unpack top_x, top_y, length, taper, dip, depth, thickness_array = slab
    current_depth = top_y
    total_thickness = sum(thickness_array)
    polyshapes = Vector{Vector{Vector{Float64}}}()
    for thickness in thickness_array
        vertices = [
            [top_x, top_y],
            [top_x + total_thickness / tand(taper), current_depth],
            [top_x + length, current_depth],
            [top_x + length + depth / tand(dip), current_depth - depth],
            [top_x + length + depth / tand(dip), current_depth - depth - thickness],
            [top_x + length, current_depth - thickness],
            [top_x + total_thickness / tand(taper), current_depth - thickness],
            [top_x, top_y],  # Closing the polygon
        ]
        push!(polyshapes, vertices)
        current_depth -= thickness
    end
    slab.polygons = polyshapes
end

# Define a function to compute vertices for Indentor2D
function compute_poly(indentor::Indentor2D)
    @unpack top_x, top_y, length, taper, taper2, thickness_array = indentor
    taper2 = taper2 === nothing ? taper : taper2
    d = top_y
    polyshapes = Vector{Vector{Vector{Float64}}}()
    for thickness in thickness_array
        vertices = [
            [top_x - d / tand(taper2), top_y + d],
            [top_x + length + d / tand(taper), top_y + d],
            [top_x + length + (d - thickness) / tand(taper), top_y + d - thickness],
            [top_x - (d - thickness) / tand(taper2), top_y + d - thickness],
            [top_x - d / tand(taper2), top_y + d],
        ]
        push!(polyshapes, vertices)
        d -= thickness
    end
    indentor.polygons = polyshapes
end

# Define a function to compute vertices for OverRidingPlate2D
function compute_poly(plate::OverRidingPlate2D)
    @unpack top_x, top_y, length, taper, dip, thickness_array = plate
    d = top_y
    polyshapes = Vector{Vector{Vector{Float64}}}()
    total_thickness = sum(thickness_array)

    for thickness in thickness_array
        vertices = [
            [top_x, top_y],
            [top_x + total_thickness / tand(taper), top_y + d],
            [top_x + length - d / tand(dip), top_y + d],
            [top_x + length + ( thickness - d) / tand(dip), top_y -  thickness + d],
            [top_x + total_thickness / tand(taper), top_y -  thickness + d],
            [top_x, top_y]
        ]
        push!(polyshapes, vertices)
        d -=  thickness
    end
    plate.polygons = polyshapes
end

# Define an outer constructor that computes the vertices
function Indentor2D(; top_x::Float64, top_y::Float64, length::Float64, taper::Float64, taper2::Union{Nothing, Float64} = nothing, thickness_array::Vector{Float64})
    i = Indentor2D(top_x, top_y, length, taper, taper2, thickness_array, nothing)
    compute_poly(i)
    return i
end

function Slab2D(; top_x::Float64, top_y::Float64, length::Float64, taper::Float64, dip::Float64, depth::Float64, thickness_array::Vector{Float64})
    s = Slab2D(top_x, top_y, length, taper, dip, depth, thickness_array, nothing)
    compute_poly(s)
    return s
end

function OverRidingPlate2D(; top_x::Float64, top_y::Float64, length::Float64, taper::Float64, dip::Float64, thickness_array::Vector{Float64})
    p = OverRidingPlate2D(top_x, top_y, length, taper, dip, thickness_array, nothing)
    compute_poly(p)
    return p
end


function Base.setproperty!(p::Shape2D, sym::Symbol, val)
    if getfield(p, sym) != val
        setfield!(p, sym, val)
        compute_polygons(p)
    end
end

function Base.setproperty!(p::Shape2D, sym::Symbol, val)
    Base.setfield!(p, sym, val)
    if sym != :polygons
        compute_poly(p)
    end
end
