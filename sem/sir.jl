using BenchmarkTools
using Plots
using StatsBase


# Check getInfectives and whether it should return the index of the infectives in ttree, or their actual identity

function getInfectives(ttree::Matrix{Float64}, t::Float64)::Vector{Int64}
    infectives = Vector{Int64}()
    for i in 1:size(ttree)[2]
        if ttree[3, i] <= t     # Check if individual has been infected
            if ttree[4, i] > t   # If so, check if they have recovered
                push!(infectives, i)
            end
        else
            break
        end
    end
    return infectives
end


struct Outbreak
    ttree::Matrix{Float64}  # Full transmission tree
    trajectory::Tuple{Vector{Float64}, Matrix{Float64}}
end

# Incorporate (immediate) offspring distribution calculation as output
function simulateOutbreak_push(R0::Float64=2.0,
                               γ::Float64=1.0,
                               nagents::Int64=1_000,
                               ninfectives::Int64=1,
                               tspan::Tuple{Float64, Float64}=(0.0, 100.0),
                               maxiter::Int64=1_000_000)::Outbreak

    # Transmission tree: [Infectee, Infector, TimeOfInfection, TimeOfRecovery, TimeOfSampling]
    ttree = fill(-1.0, (5, nagents))  # Initialize a transmission tree
    ttree[1, 1:ninfectives] .= collect(1:ninfectives)
    ttree[2,1:ninfectives] .= fill(0.0, ninfectives)
    ttree[3,1:ninfectives] .= fill(tspan[1], ninfectives)
    ttree[4,1:ninfectives] .= fill(Inf, ninfectives)
    cumincidence = ninfectives
    t = tspan[1]
    infectives = getInfectives(ttree, t)

    nagents = convert(Float64, nagents)
    ninfectives = convert(Float64, ninfectives)

    S = nagents - ninfectives
    I = ninfectives
    R = 0.0

    tout = [t]
    Sout = [S]
    Iout = [I]
    Rout = [R]

    iter = 0

    while (length(infectives) > 0) && (iter <= maxiter) && (t <= tspan[2])
        iter += 1        

        infection_rate = R0 * γ * S * I / nagents
        recovery_rate = γ * I

        r = rand()
        t -= log(r) / (infection_rate + recovery_rate)
        if r <= infection_rate / (infection_rate + recovery_rate)
            infector = sample(infectives)
            cumincidence += 1
            ttree[:, cumincidence] .= [convert(Float64, cumincidence), infector, t, Inf, -1.0]
            S -= 1.0
            I += 1.0
        else
            recovered = sample(infectives)
            ttree[4, Int(recovered)] = t
            I -= 1.0
            R += 1.0
        end
        push!(tout, t)
        push!(Sout, S)
        push!(Iout, I)
        push!(Rout, R)
        infectives = getInfectives(ttree, t)

        end
        return Outbreak(ttree, (tout, hcat(Sout, Iout,Rout)))
end






function simulateOutbreak(R0::Float64=2.0,
                          γ::Float64=1.0,
                          nagents::Int64=1_000,
                          ninfectives::Int64=1,
                          tspan::Tuple{Float64, Float64}=(0.0, 100.0),
                          maxiter::Int64=1_000_000)::Matrix{Float64}
    # Transmission tree: [Infectee, Infector, TimeOfInfection, TimeOfRecovery, TimeOfSampling]
    ttree = fill(-1.0, (5, nagents))  # Initialize a transmission tree
    ttree[1, 1:ninfectives] .= collect(1:ninfectives)
    ttree[2,1:ninfectives] .= fill(0.0, ninfectives)
    ttree[3,1:ninfectives] .= fill(tspan[1], ninfectives)
    cumincidence = ninfectives
    infectives = getInfectives(ttree)

    iter = 0
    t = tspan[1]

    while (length(infectives) > 0) && (iter <= maxiter) && (t <= tspan[2])
        iter += 1        

        infection_rate = R0 * γ * (nagents - cumincidence) * length(infectives) / nagents
        recovery_rate = γ * length(infectives)

        r = rand()
        t -= log(r) / (infection_rate + recovery_rate)
        if r <= infection_rate / (infection_rate + recovery_rate)
            infector = sample(infectives)
            cumincidence += 1
            ttree[:, cumincidence] .= [cumincidence, infector, t, -1.0, -1.0]
        else
            recovered = sample(infectives)
            ttree[4, Int(recovered)] = t
        end

        infectives = getInfectives(ttree)

    end
    return ttree
end


function outbreakSize(ttree::Matrix{Float64})::Int
    os = 0
    for col in eachcol(ttree)
        if col[1] > 0.0
            os += 1
        else
            break
        end
    end
    return os
end


function sampleSize(ttree::Matrix{Float64})::Int
    ss = 0
    for col in eachcol(ttree)
        if col[5] >= 0.0
            ss += 1
        else
           break 
        end
    end
    return ss
end


function getSampled(ttree::Matrix{Float64})::Vector{Float64}
    sampled = Vector{Float64}()
    for col in eachcol(ttree)
        if col[5] >= 0.0
            push!(sampled, col[1])
        end
    end
    return sampled
end

# This function may be pretty inefficient
function toDict(ttree::Matrix{Float64})::Dict{Float64, Tuple{Float64, Float64, Float64, Float64}}
    dct = Dict{Float64, Tuple{Float64, Float64, Float64, Float64}}()
    for col in eachcol(ttree)
        if col[1] >= 0.0
            dct[col[1]] = (col[2], col[3], col[4], col[5])
        else
            break
        end
    end
    return dct
end


function sampleOutbreak(ttree::Matrix{Float64},
                        T::Float64,
                        π::Float64)::Matrix{Float64}
    ttree_ = copy(ttree)
    sampled = Vector{Float64}()
    infectives = getInfectives(ttree_, T)
    for i in infectives
        if rand() < π
            push!(sampled, i)
        end
    end
    # Set sampling time of sampled infectives
    ttree_[5, Int.(sampled)] .= T
    # Now traverse the transmission tree folowing the transmission chains from infectee to infector
    # for each sampled infective
    leaves = copy(sampled)
    while (length(leaves) > 0)
        infector = ttree[2, Int(pop!(leaves))]
        if !(infector in sampled) && infector > 0.0
            push!(sampled, infector)
            push!(leaves, infector)
        end
    end
    return ttree_[:, Int.(sampled)]
end



infectees = @view s[1, :]
infectors = @view s[2, :]

infecteds = sort(unique(vcat(infectess, infectors)), rev=true)

directoffspring = Dict(zip(infecteds, [Vector{Float64}() for _ in infecteds]))

for col in eachcol(s)
    push!(directoffspring[col[2]], col[1])
end

totprogeny = Dict(zip(infecteds, [0 for _ in infecteds]))

for i in infecteds
    if length(directoffspring[i]) == 0
        totprogeny[i] += 1
    else
        for off in directoffspring[i]
            totprogeny[i] += totprogeny[off]
        end
    end
end

infectives = getInfectives(s, 2.0)




function offspringDist(ttree::Matrix{Float64})::Dict{Float64, Int64}
    sampled = getSampled(ttree)


end






function newick(ttree::Matrix{Float64})
    
    sampled = getSampled(ttree)
    if length(sampled) <= 1
        return nothing
    end
    leaves = copy(sampled)
    stree = ttree[:, sortperm(ttree[3, :])]

    while length(leaves) > 0
        


    

    # Convert array to dictionary
    tree = Dict{Float64, Tuple{Float64, Tuple{Float64, Float64, Float64}}}()
    for col in eachcol(ttree)
        tree[col[1]] = (col[2], (col[3], col[4], col[5]))
    end

    intnode = sampled + 1
    child = pop!(sampled)
    parent = tree[child]
    if parent in sampled
        



# This function is very inefficient!
function epidemicTrajectory(ttree::Matrix{Float64})::Matrix{Float64}
    N = size(ttree)[2]
    t0 = ttree[2,1]
    tinfections = filter(x -> x >= t0, ttree[3, :]) # These values should already be sorted
    trecoveries = filter(x -> x >= t0, ttree[4, :])
    tevents = sort(vcat(tinfections, trecoveries))
    
    S = fill(-1.0, length(tevents))
    I = fill(-1.0, length(tevents))
    R = fill(-1.0, length(tevents))

    S[1] = N - 1.0
    I[1] = 1.0
    R[1] = 0.0

    for idx in 2:length(tevents)
        if tevents[idx] in tinfections
            S[idx] = S[idx-1] - 1.0
            I[idx] = I[idx-1] + 1.0
            R[idx] = R[idx-1]
        else
            S[idx] = S[idx-1]
            I[idx] = I[idx-1] - 1.0
            R[idx] = R[idx-1] + 1.0
        end        
    end
    return(hcat(tevents, S, I, R))
end



ttree = simulateOutbreak();
out = epidemicTrajectory(ttree);
plot(out[])