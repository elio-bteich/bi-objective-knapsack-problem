struct dataKP
    p::Vector{Float64} # profits objectif
    w::Vector{Float64} # poids
    Omega::Int64 # capacité
    nbItems::Int64 # nombre d'objets
end

struct data2KP
    p1::Vector{Int64} # profits objectif 1
    p2::Vector{Int64} # profits objectif 2
    w::Vector{Int64} # poids dimension 1
    Omega::Int64 # capacité dimension 1
    nbItems::Int64 # nombre d'objets
end

mutable struct CriterionVector2KP
    z1::Int64
    z2::Int64
    pred::Union{CriterionVector2KP, Nothing}
end

const State2KP = Vector{CriterionVector2KP}

mutable struct solutionKP
    x::Vector{Int8}
    z1::Int64
    z2::Int64
end

function parseKP(filename::String)
	f::IOStream = open(filename,"r")

    # Première ligne : taille du problème (nombre d'objets)
    s::String = readline(f)
    tab::Vector{Int64} = parse.(Int64,split(s," ",keepempty = false))
    @inbounds n::Int64 = tab[1]

	# Deuxième ligne : capacité du sac à dos
    s = readline(f)
	tab = parse.(Int64,split(s," ",keepempty = false))
	@inbounds Omega::Int64 = tab[1]

	# Troisième ligne : coefficients des coefficients de la première fonction objectif
    s = readline(f)
	p1::Vector{Int64} = parse.(Int64,split(s," ",keepempty = false))

    # Quatrième ligne : coefficients des coefficients de la seconde fonction objectif
    s = readline(f)
	p2::Vector{Int64} = parse.(Int64,split(s," ",keepempty = false))

	# Cinquième ligne : poids des objets
    s = readline(f)
	w::Vector{Int64} = parse.(Int64,split(s," ",keepempty = false))

    # End
    close(f)

    return data2KP(p1,p2,w,Omega,n)
end

