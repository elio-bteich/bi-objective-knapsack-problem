using Plots
include("./struct_parse.jl")

# Plot a list of vectors on a plane to visualize better
function plot_vectors(vectors::State2KP)
    x::Vector{Int64} = []
    y::Vector{Int64} = []
    
    for vector in vectors
        push!(x, vector.z1)
        push!(y, vector.z2)
    end

    scatter(x, y)
    
end

# Plot a list of vectors on a plane to visualize better
function plot_vectors(vectors::Vector{Tuple{Int64, Int64}})
    x::Vector{Int64} = []
    y::Vector{Int64} = []
    
    for vector in vectors
        push!(x, vector[1])
        push!(y, vector[2])
    end

    scatter(x, y)
    
end