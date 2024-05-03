using Plots
include("./struct_parse.jl")

# Plot a state (criterion vectors) on a plane. (Only compatible with new version of prog-dyn)
function plot_vectors(vectors::State2KP, bound_points::Vector{Tuple{Float64, Float64}}, vectors_desc::String, bound_desc::String)
    x::Vector{Int64} = []
    y::Vector{Int64} = []

    # Extract x and y coordinates from the original vectors
    for vector in vectors
        push!(x, vector.z1)
        push!(y, vector.z2)
    end

    # Create the initial scatter plot
    p = scatter(x, y, label=vectors_desc, color=:blue, markersize=4)

    # Extract x and y coordinates for bound_points
    # Sort the bound_points by the first coordinate of each tuple
    sorted_bound_points = sort(bound_points, by=first)

    bx::Vector{Float64} = Float64[]
    by::Vector{Float64} = Float64[]
    
    for point in sorted_bound_points
        push!(bx, point[1])
        push!(by, point[2])
    end

    plot!(bx, by, color=:red, label=bound_desc, lw=4)


    return p
end



# Plot a state (criterion vectors) on a plane. (Only compatible with old version of prog-dyn)
function plot_vectors(vectors::Vector{Tuple{Int64, Int64}})
    x::Vector{Int64} = []
    y::Vector{Int64} = []
    
    for vector in vectors
        push!(x, vector[1])
        push!(y, vector[2])
    end

    scatter(x, y)
    
end