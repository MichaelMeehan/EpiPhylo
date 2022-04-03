using BenchmarkTools
using Plots
using StatsBase


struct Outbreak2
    tree::Matrix{Float64}  # Transmission tree: array of [Infectee id; Infector id; Infection time; Recovery time; Sampling time]
    off::Dict{Float64, Vector{Tuple{Float64, Float64}}}     # Offspring distribution: dictionary of who infected whom and time of infection
    anc::Dict{Float64, Tuple{Float64, Float64}}     # Ancestry distribution: dictionary of ancestry of each infected individual
    traj::Tuple{Vector{Float64}, Matrix{Float64}}   # Epidemic trajectory (S(t), I(t), R(t))
end


function getInfectives(tree::Matrix{Float64}, t::Float64)::Vector{Float64}
    infectives = Vector{Float64}()
    for col in eachcol(tree)
        if col[1] >= 0.0 && col[3] <= t
            if col[4] > t
                push!(infectives, col[1])
            end
        else
            break
        end
    end
    return infectives
end


function simulateOutbreak2(R0::Float64=2.0,
                          γ::Float64=1.0,
                          nagents::Int64=1_000,
                          ninfectives::Int64=1,
                          tspan::Tuple{Float64, Float64}=(0.0, 100.0),
                          maxiter::Int64=1_000_000)::Outbreak2

    # Transmission tree: [Infectee, Infector, TimeOfInfection, TimeOfRecovery, TimeOfSampling]
    tree = fill(-1.0, (5, nagents))  # Initialize a transmission tree
    infectives = collect(1.0:ninfectives)
    tree[1, 1:ninfectives] .= infectives
    tree[2,1:ninfectives] .= fill(0.0, ninfectives)
    tree[3,1:ninfectives] .= fill(tspan[1], ninfectives)
    tree[4,1:ninfectives] .= fill(Inf, ninfectives)

    cumincint = ninfectives
    t = tspan[1]

    off = Dict(zip(infectives, [Vector{Tuple{Float64, Float64}}() for _ in infectives]))
    off[0.0] = [(i, t) for i in infectives]
    
    anc = Dict(zip(infectives, [(0.0, t) for _ in infectives]))

    nagents = convert(Float64, nagents)
    ninfectives = convert(Float64, ninfectives)
    cumincfloat = ninfectives

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
            cumincint += 1
            cumincfloat += 1
            tree[:, cumincint] .= [cumincfloat, infector, t, Inf, -1.0]

            if infector in keys(off)
                push!(off[infector], (cumincfloat, t))
            else
                off[infector] = [(cumincfloat, t)]
            end

            anc[cumincfloat] = (infector, t)

            S -= 1.0
            I += 1.0
        else
            recovered = sample(infectives)
            tree[4, Int(recovered)] = t
            I -= 1.0
            R += 1.0
        end
        push!(tout, t)
        push!(Sout, S)
        push!(Iout, I)
        push!(Rout, R)
        infectives = getInfectives(tree, t)

        end
        return Outbreak2(tree, off, anc, (tout, hcat(Sout, Iout,Rout)))
end


function sampleOutbreak(out::Outbreak,
                        t::Float64,
                        π::Float64)::Outbreak
    stree = copy(out.tree)
    anc = Dict{Float64, Tuple{Float64, Float64}}()
    retained = Vector{Float64}()
    infectives = getInfectives(out.tree, t)
    for i in infectives
        if rand() < π
            stree[5, Int(i)] = t
            ancestor = i
            while ancestor > 0.0
                if !(ancestor in retained)
                    anc[ancestor] = out.anc[ancestor]
                    push!(retained, ancestor)
                    ancestor = out.anc[ancestor][1]
                else
                    break
                end
            end
        end
    end
    
    off = Dict{Float64, Vector{Tuple{Float64, Float64}}}()
    for ret in vcat(0.0, retained)
        off[ret] = Vector{Tuple{Float64, Float64}}()
        if ret in keys(out.off)
            for itpair in out.off[ret]
                if itpair[1] in retained
                    push!(off[ret], itpair)
                end
            end
        end
    end
    return Outbreak(stree[:, sort(Int.(retained))], off, anc, out.traj)
end


function getSampled(out::Outbreak)::Vector{Float64}
    sampled = Vector{Float64}()
    for col in eachcol(out.tree)
        if col[5] >= 0.0
            push!(sampled, col[1])
        end
    end
    return sampled
end


function getClusters(out::Outbreak)::Tuple{Vector{Float64}, Vector{Dict{Float64, Vector{Float64}}}}
    # Collect list of sampled individuals
    sampled = getSampled(out)
    # Generate cluster-size distribution backward in time
    clusters = Dict(zip(sampled, [[s] for s in sampled]))
    t = out.tree[5, end]
    tout = [t]
    cout = [clusters]
    for col in size(out.tree)[2]:-1:1
        infectee = out.tree[1, col]
        infector = out.tree[2, col]
        t = out.tree[3, col]
        leaves = keys(clusters)
        # Check if two nodes are leaf nodes
        if infector in leaves
            # Merge clusters
            clusters[infector] = vcat(clusters[infectee], clusters[infector])
        else
            # Update ancestor
            clusters[infector] = clusters[infectee]
        end
        delete!(clusters, infectee)
        push!(tout, t)
        push!(cout, copy(clusters))
    end
    cout[1] = Dict(zip(sampled, [[s] for s in sampled]))
    return tout, cout
end


function getCSDVec(clusters::Vector{Dict{Float64, Vector{Float64}}})::Vector{Vector{Int64}}
    csdvec = Vector{Vector{Int}}() # Inititate an empty vector to store cluster size distributions
    for ec in clusters   # Loop over the sequence of equivalence classes
        dist = Vector{Int64}() 
        for c in values(ec) # Loop over each cluster in the equivalence class
            push!(dist, length(c))
        end
        push!(csdvec, dist)
    end
    return csdvec
end


function getCSD(csdvec::Vector{Vector{Int64}}, tinf::Vector{Float64}, t::Float64)::Vector{Int64}
    if t >= tinf[1]
        return csdvec[1]
    else
        for idx in eachindex(tinf)
            if tinf[idx] <= t
                return csdvec[idx-1]
            end
        end
    end
    return csdvec[end]
end

# T = [7.96, 9.77, 11.22, 16.06]
# t = [[5.27, 0.96, 2.68, 0.1, 6.13, 6.99, 3.54, 7.86, 4.41, 1.82], [4.35, 0.1, 9.67, 5.42, 3.29, 1.16, 8.61, 6.48, 2.22, 7.54], [4.99, 0.1, 6.22, 1.32, 8.67, 3.77, 9.89, 7.44, 11.12, 2.54], [3.62, 5.38, 0.1, 7.15, 1.86, 8.91, 10.67, 12.44, 14.20, 15.96]]

out = simulateOutbreak(1.0/0.3, 0.3, 100_000, 1_000);

# plot(out.traj[1][2:end], out.traj[2][2:end,2]/100_000, xaxis=:log10, xlims=(0.1, 20), grid=true, minorticks=10)

s = sampleOutbreak(out, 7.96, 1.0);

tinf, clusters = getClusters(s);

csdvec = getCSDVec(clusters);

csds = map(t -> getCSD(csdvec, tinf, t), sort([5.27, 0.96, 2.68, 0.1, 6.13, 6.99, 3.54, 7.86, 4.41, 1.82]))

csd = getCSD(csdvec, tinf, 0.1);

mean(csd)

tvec = out.traj[1]
I = out.traj[2][:,2]

A = [length(getInfectives(s.tree, t)) / 100_000 for t in 0.1:0.01:7.96]

x = length(getInfectives(s.tree, 7.96)) ./ A