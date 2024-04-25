# Author: Elio Bteich
# Date: 17 April 2024

include("./util.jl")
include("./struct_parse.jl")

# Return true if v1 dominates v2, false otherwise
function is_dominant(v1::CriterionVector2KP, v2::CriterionVector2KP)

    return v1.z1 >= v2.z1 && v1.z2 >= v2.z2 

end

# Find the non-dominated criterion vectors between 2 states of non-dominated vectors
function non_dominated_vectors(S1::State2KP, S2::State2KP, elemIndex::Int64, d::data2KP)

    nd::State2KP = []

    p_s1 = 1
    p_s2 = 1
    p_nd = 0

    # Since the nd vector is empty, we have treated the first iteration apart from the others to prevent checking if the nd vector is empty everytime we want to insert a new non-dominant vector in it
    if S1[p_s1].z1 > S2[p_s2].z1 + d.p1[elemIndex]    

        push!(nd, S1[p_s1])
        p_s1 += 1
        
    elseif S1[p_s1].z1 < S2[p_s2].z1 + d.p1[elemIndex]
        
        push!(nd, CriterionVector2KP(S2[p_s2].z1 + d.p1[elemIndex], S2[p_s2].z2 + d.p2[elemIndex], S2[p_s2]))
        p_s2 += 1

    else
        
        if S1[p_s1].z2 >= S2[p_s2].z2 + d.p2[elemIndex]
            
            push!(nd, CriterionVector2KP(S1[p_s1].z1, S1[p_s1].z2, S1[p_s1]))

        else
            
            push!(nd, CriterionVector2KP(S2[p_s2].z1 + d.p1[elemIndex], S2[p_s2].z2 + d.p2[elemIndex], S2[p_s2]))            

        end
        
        p_s1 += 1
        p_s2 += 1
                    
    end

    p_nd += 1
    
    while p_s1 <= size(S1, 1) && p_s2 <= size(S2, 1)

        if S1[p_s1].z1 > S2[p_s2].z1 + d.p1[elemIndex]
        
            if nd[p_nd].z2 < S1[p_s1].z2

                push!(nd, CriterionVector2KP(S1[p_s1].z1, S1[p_s1].z2, S1[p_s1]))
                p_nd += 1

            end
            
            p_s1 += 1
            
        elseif S1[p_s1].z1 < S2[p_s2].z1 + d.p1[elemIndex]
            
            if nd[p_nd].z2 < S2[p_s2].z2 + d.p2[elemIndex]
                
                push!(nd, CriterionVector2KP(S2[p_s2].z1 + d.p1[elemIndex], S2[p_s2].z2 + d.p2[elemIndex], S2[p_s2]))
                p_nd += 1

            end
            
            p_s2 += 1
        
        else
            
            if S1[p_s1].z2 >= S2[p_s2].z2 + d.p2[elemIndex]
                
                if nd[p_nd].z2 < S1[p_s1].z2

                    push!(nd, CriterionVector2KP(S1[p_s1].z1, S1[p_s1].z2, S1[p_s1]))
                    p_nd += 1

                end
            
            else
                
                if nd[p_nd].z2 < S2[p_s2].z2 + d.p2[elemIndex]
            
                    push!(nd, CriterionVector2KP(S2[p_s2].z1 + d.p1[elemIndex], S2[p_s2].z2 + d.p2[elemIndex], S2[p_s2]))
                    p_nd += 1

                end
            
            end
            
            p_s1 += 1
            p_s2 += 1
                        
        end

    end

    # S1 still has vectors that we haven't iterated on
    while p_s1 <= size(S1, 1)
    
        if nd[p_nd].z2 < S1[p_s1].z2

            push!(nd, CriterionVector2KP(S1[p_s1].z1, S1[p_s1].z2, S1[p_s1]))
            p_nd += 1
        end

        p_s1 += 1
        
    end
   
    # S2 still has vectors that we haven't iterated on
    while p_s2 <= size(S2, 1)

        if nd[p_nd].z2 < S2[p_s2].z2 + d.p2[elemIndex]

            push!(nd, CriterionVector2KP(S2[p_s2].z1 + d.p1[elemIndex], S2[p_s2].z2 + d.p2[elemIndex], S2[p_s2]))
            p_nd += 1
        end

        p_s2 += 1
        
    end

    return nd

end


# Calculate the states matrix of the multi-objective knapsack problem, for line i and column k each state contains the non-dominated criterion vectors of the problem with the first k elements and the capacity i
function calculate_matrix(d::data2KP)

    if d.nbItems <= 0
        return Matrix(undef, 0, 0)
    end

    R::Matrix{State2KP} = Matrix{State2KP}(undef, d.Omega + 1, d.nbItems)

    # column 1    
    for i in 1:d.w[1]

        R[i, 1] = [CriterionVector2KP(0, 0, nothing)]
   
    end

    for i in (d.w[1]+1):(d.Omega+1)
    
        R[i, 1] = [CriterionVector2KP(d.p1[1], d.p2[1], nothing)]

    end

    # other columns
    for j in 2:d.nbItems
        
        for i in 1:d.w[j]

            n = length(R[i, j-1])

            R[i, j] = State2KP(undef, n)

            for t in 1:n
            
                R[i, j][t] = CriterionVector2KP(R[i, j-1][t].z1, R[i, j-1][t].z2, R[i, j-1][t])

            end
        
        end

        for i in (d.w[j]+1):(d.Omega+1)

            R[i, j] = non_dominated_vectors(R[i, j-1], R[i-d.w[j], j-1], j, d)
        
        end
    
    end

    return R

end

# Extracts the non-dominated criterion vectors from the last cell of the matrix, which corresponds to the final solution space.
function find_nd_criterion_vectors(R::Matrix{State2KP})

    return R[end, end]

end

function solve_prog_dyn(d::data2KP)

    R = calculate_matrix(d)
    return find_nd_criterion_vectors(R)

end

function solve_v3()

    d = parseKP("Instances KP/2KP100-10.dat")
    @time nd::State2KP = solve_prog_dyn(d)
    println(length(nd))
    plot_vectors(nd)

end


