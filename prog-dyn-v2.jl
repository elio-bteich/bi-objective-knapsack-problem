# Author: Elio Bteich
# Date: 17 April 2024

include("./util.jl")
include("./struct_parse.jl")

# Return true if v1 dominates v2, false otherwise
function is_dominant(v1::Tuple{Int64, Int64}, v2::Tuple{Int64, Int64})

    return v1[1] >= v2[1] && v1[2] >= v2[2] 

end

# Find the non-dominated criterion vectors between 2 states of non-dominated vectors
function non_dominated_vectors(S1::Vector{Tuple{Int64, Int64}}, S2::Vector{Tuple{Int64, Int64}})

    nd::Vector{Tuple{Int64, Int64}} = []

    p_s1 = 1
    p_s2 = 1
    p_nd = 0

    # Since the nd vector is empty, we have treated the first iteration apart from the others to prevent checking if the nd vector is empty everytime we want to insert a new non-dominant vector in it
    if S1[p_s1][1] > S2[p_s2][1]        

        push!(nd, S1[p_s1])
        p_s1 += 1
        p_nd += 1
        
    elseif S1[p_s1][1] < S2[p_s2][1]
        
        push!(nd, S2[p_s1])
        p_s2 += 1
        p_nd += 1

    else
        
        if S1[p_s1][2] >= S2[p_s2][2]
            
            push!(nd, S1[p_s1])
            p_nd += 1

        else
            
            push!(nd, S2[p_s2])
            p_nd += 1

        end
        
        p_s1 += 1
        p_s2 += 1
                    
    end

    
    while p_s1 <= size(S1, 1) && p_s2 <= size(S2, 1)

        if S1[p_s1][1] > S2[p_s2][1]
        
            if nd[p_nd][2] < S1[p_s1][2]

                push!(nd, S1[p_s1])
                p_nd += 1

            end
            
            p_s1 += 1
            
        elseif S1[p_s1][1] < S2[p_s2][1]
            
            if nd[p_nd][2] < S2[p_s2][2]
                
                push!(nd, S2[p_s2])
                p_nd += 1

            end
            
            p_s2 += 1
        
        else
            
            if S1[p_s1][2] >= S2[p_s2][2]
                
                if nd[p_nd][2] < S1[p_s1][2]

                    push!(nd, S1[p_s1])
                    p_nd += 1

                end
            
            else
                
                if nd[p_nd][2] < S2[p_s2][2]
            
                    push!(nd, S2[p_s2])
                    p_nd += 1

                end
            
            end
            
            p_s1 += 1
            p_s2 += 1
                        
        end

    end

    # S1 still has vectors that we haven't iterated on
    while p_s1 <= size(S1, 1)
    
        if nd[p_nd][2] < S1[p_s1][2]

            push!(nd, S1[p_s1])
            p_nd += 1
        end

        p_s1 += 1
        
    end
   
    # S2 still has vectors that we haven't iterated on
    while p_s2 <= size(S2, 1)

        if nd[p_nd][2] < S2[p_s2][2]

            push!(nd, S2[p_s2])
            p_nd += 1
        end

        p_s2 += 1
        
    end

    return nd

end

# Returns the criterion vectors after the addition of an element
function vectors_with_element(vectors::Vector{Tuple{Int64, Int64}}, elemIndex::Int64, d::data2KP)

    vectorAdded::Tuple{Int64, Int64} = (d.p1[elemIndex], d.p2[elemIndex]) 

    n = size(vectors, 1)

    vectorsWithElem::Vector{Tuple{Int64, Int64}} = Vector{Tuple{Int64, Int64}}(undef, n)

    for i in 1:n
        vectorsWithElem[i] = (vectors[i][1] + vectorAdded[1], vectors[i][2] + vectorAdded[2])
    end

    return vectorsWithElem

end

# Calculate the states matrix of the multi-objective knapsack problem, for line i and column k each state contains the non-dominated criterion vectors of the problem with the first k elements and the capacity i
function calculate_matrix(d::data2KP)

    if d.nbItems <= 0
        return Matrix(undef, 0, 0)
    end

    R::Matrix{Vector{Tuple{Int64, Int64}}} = Matrix{Vector{Tuple{Int64, Int64}}}(undef, d.Omega + 1, d.nbItems)

    # column 1    
    for i in 1:d.w[1]

        R[i, 1] = [(0,0)]
   
    end

    for i in (d.w[1]+1):(d.Omega+1)
    
        R[i, 1] = [(d.p1[1], d.p2[1])]

    end

    # other columns
    for j in 2:d.nbItems
        
        for i in 1:d.w[j]
        
            R[i, j] = R[i, j-1]

        end

        for i in (d.w[j]+1):(d.Omega+1)

            R[i, j] = non_dominated_vectors(vectors_with_element(R[i-d.w[j], j-1], j, d), R[i, j-1])
        
        end
    
    end

    return R

end

# Find the non-dominant criterion vectors from the states matrix of the dynamic programming algorithm
function find_nd_criterion_vectors(R::Matrix{Vector{Tuple{Int64, Int64}}})
    return R[end, end]
end

function solve_prog_dyn(d::data2KP)

    R = calculate_matrix(d)
    return find_nd_criterion_vectors(R)

end

function solve_v2()

    d = parseKP("Instances KP/2KP100-10.dat")
    @time nd::Vector{Tuple{Int64, Int64}} = solve_prog_dyn(d)
    println(length(nd))
    plot_vectors(nd)

end



