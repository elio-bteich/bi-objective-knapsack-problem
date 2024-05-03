# Author: Elio Bteich
# Date: 30 April 2024

include("./util.jl")
include("./struct_parse.jl")

# Return the permutation that sorts elements in descending order of the profit/weight ratio according to the weighted sum
function get_weighted_sum_sorting_permutation(lambda1::Float64, lambda2::Float64, data::data2KP)

    if (lambda1 == 0 && lambda2 != 0)

        return sortperm([1000*data.p2[i] / data.w[i] + data.p1[i] / data.w[i] for i in 1:data.nbItems], rev=true)

    elseif (lambda1 != 0 && lambda2 == 0)

        return sortperm([1000*data.p1[i] / data.w[i] + data.p2[i] / data.w[i] for i in 1:data.nbItems], rev=true)

    end

    return sortperm([lambda1 * data.p1[i] / data.w[i] + lambda2 * data.p2[i] / data.w[i] for i in 1:data.nbItems], rev=true)

end

# Calculate the continuous relaxation of the problem combining linearly the 2 objective functions of the bi-objective problem
function continuous_relaxation(lambda1::Float64, lambda2::Float64, data::data2KP)

    perm::Vector{Int64} = get_weighted_sum_sorting_permutation(lambda1, lambda2, data)
    
    sorted_p1 = data.p1[perm]
    sorted_p2 = data.p2[perm]
    sorted_w = data.w[perm]

    residual_capacity::Float64 = float(data.Omega)
    i::Int64 = 1
    z1::Float64 = 0
    z2::Float64 = 0

    while (residual_capacity > 0 && i <= data.nbItems)

        if (sorted_w[i] <= residual_capacity)
            residual_capacity -= sorted_w[i]
            z1 += sorted_p1[i]
            z2 += sorted_p2[i]
        else
            residual_ratio::Float64 = residual_capacity / sorted_w[i]
            z1 += residual_ratio * sorted_p1[i]
            z2 += residual_ratio * sorted_p2[i]
            residual_capacity = 0
        end

        i += 1

    end

    return z1, z2

end

# Caclulate the continuous relaxation bound of a bi-objective problem with the weighted sum method
function calculate_bound(data::data2KP)

    y1 = continuous_relaxation(Float64(0.0), Float64(1.0), data)
    y2 = continuous_relaxation(Float64(1.0), Float64(0.0), data)

    bound::Vector{Tuple{Float64, Float64}} = [y1, y2]

    calculate_bound_rec(y1, y2, bound, data)

    return bound
end

# Caclulate the continuous relaxation bound of a bi-objective problem with the weighted sum method
function calculate_bound_rec(y1::Tuple{Float64, Float64}, y2::Tuple{Float64, Float64}, bound::Vector{Tuple{Float64, Float64}}, data::data2KP)

    lambda1::Float64 = y1[2] - y2[2]
    lambda2::Float64 = y2[1] - y1[1]

    y::Tuple{Float64, Float64} = continuous_relaxation(lambda1, lambda2, data)

    push!(bound, y)

    if lambda1*(y[1] - y1[1]) > lambda2*(y1[2] - y[2])

        calculate_bound_rec(y1, y, bound, data)
        calculate_bound_rec(y, y2, bound, data)

    end

end
