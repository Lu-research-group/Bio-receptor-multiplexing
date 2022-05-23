######################################################################
#                                                                    #
# Overdamped Markovian Langevin dynamics for ligand binding events   #
#     Author: Sa Hoon Min, Asawari Pagare, Zhiyue Lu                 #
#     Date: September 13, 2021 (v0.3.1)                              #
#                                                                    #
#####################################################################

#=
Version 0.5.1
    - new run_binding function with ls and nth arguments.
        It doesn't provide output instead probability.

Version 0.5.0
    - Without WCA potential


Version 0.4.3
    - Optimization
        to improve performance further,
        1. pos_lig -> pos_x, pos_y, pos_z
        2. force -> force_x, force_y, force_z
        3. in update_force! for k and j loop -> kloop and jloop vectors

Version 0.4.1
    - DEBUG!!! 
        update_force! was revised. Previsouly, it has wrong initialization.
    - Optimization
        1. max(1,2,3) -> if 1 <= cutoff && 2 <= cutoff && 3 <= cutoff
        2. @inbounds for loops
    - Energy minimization added, steepest_descent! function
        System may have very high energy state because of overlap.
        So, before the main calculation, the random initial coordinates
        need to be relaxed.
    - calc_closest_lig
        To avoid unphysical contact between ligands during the binding event,
        now the receptor binds with the closest ligand not a random ligand, 
        when there are more than two ligands in the cutoff of receptor.

Version 0.4.0
    - Clean up

Version 0.3.11
    - new functions: trj2moreseg, callcLff (previous one is now old_calcLLF)
                     freq2Sn (removed)

Version 0.3.10
    - Updated sampling function to sample from frequency array
Version 0.3.9
    - Minor corrections

Version 0.3.8
    - Corrected input parameters
    - Removed buildHessian function; it is included in a seperate file
Version 0.3.7
    - MLE functions added
Version 0.3.6
    - In the run_binding function for making trajectories, the position 
      paramter of ligands was set to the position of rectpor, which is wrong.
      I changed "pos_rec" to "pos_lig".
    - Nr is not used here anymore, because Nr is 1. However, conventaionally,
      I left the Nr for the input paramter.

Version 0.3.5
    - Fixed an error related to idx_lig.
      Previous, idx_lig was a local variable, so its value was not updated
      at all. Internally, the value of idx_lig was always 0. So the ligand
      cannot maintain its binding state.
      To solve this local variable problem, idx_lig is now a vector whose length
      is 1. i.e., idx_lig = [0], so we can access the index of ligand through
      idx_lig[1].

Version 0.3.4
    - analysis with vector output

Version 0.3.2
    - minor polishing (removing unnecessary input parameters,
                       adding vx parameter into run_binding for writing trj)

Version 0.3.1
    - d_rec and m were removed
    - output_bit and cumsum_bit were removed
    - pos_rec set to box center initially
Version 0.3.0
    - Code changed from overdamped LD to Markovian overdamped LD
    - attractive potential between ligand and receptor was removed
    - Temp was made an input parameter
    - New input parameters were added: cutoff(cutoff r to determine binding),
      Eon(Energy of binding), Eoff(Energy of unbinding), Eb(Energy of barrier),
      alpha(to be used in drift calculation)
    - New variable idx_lig added
        -0 if receptor is unbound;
        -index of ligand bound to the receptor
    - New core function update_idx_output! added:
        -This takes care of binding, unbinding events and gives output
            - P_on = dt*(No. of ligands within cutoff)*Exp(-(Eb-Eoff)/T)
            - P_off = dt*Exp(-(Eb-Eon-gamma*vx*alpha)/T)

    - Core function update_force_output! changed to update_force!
        - Calculates force between ligands
    - Core function update_pos_rec! changed to update_pos_rec_idx!
        -if bound receptor then pos_lig of the bound receptor is same as pos_lig

Version 0.2.4:
    - cumsum_bit was removed (an output of the run_binding function)

    - waiting time calcution was added

      e.g.,
      b = calc_wait(a.output_bit)

      b.wt_b_u        : raw data of waiting time from binding to unbinding
      b.wt_u_b        : raw data of waiting time from unbinding to binding
      b.wtime_output  : [1st moment of wt_b_u, 2nd moment of wt_b_u,
                         3rd moment of wt_b_u, 1st moment of wt_u_b,
                         2nd moment of wt_u_b, 3rd moment of wt_u_b]

Version 0.2.2:
    - triple autocorrelation function was added.

      make_acf3

             ⟨f(x)f(x+τ)f(x+2τ)⟩
      h(τ) = ------------------
                ⟨f(x)f(x)⟩

Version 0.2.1:
    - Autocorrelation function (ACF) was revised.
      The previous definition of ACF in the autocor (StatsBase module) is that,

              Σ((f(x) - ⟨f(x)⟩) * (f(x+τ) - ⟨f(x)⟩))
      H(τ) = --------------------------------------
              Σ((f(x) - ⟨f(x)⟩) * (f(x) - ⟨f(x)⟩))

      If we don't use the subtraction of mean,
      e.g., b = autocor(x,lags, demean=false)
      then the definition is as follows,

              Σ(f(x)f(x+τ))
      HH(τ) = --------------
               Σ(f(x)f(x))


      The new definition of ACF (in the make_acf_dot function) is
              ⟨f(x)f(x+τ)⟩
      h(τ) = -------------
               ⟨f(x)f(x)⟩

    - make_acf_dot and make_acf_forloop are equivalent, but
      make_acf_dot is 2-3 faster than make_acf_forloop.
      Please use make_acf_dot.
      e.g., b = make_acf_dot(a.output_bit, lags)

    - Now, run_binding function will run equilibrium step (1 million steps)
      first, then will run "Nt" steps for production run.
      So, there is no "output_cut" data anymore.

    - If one doesn't need to make output files of trajectory and output.txt,
      just execute the run_binding function without the freq and o_name.

      e.g., run_binding(Nl,Nr,box_x,box_y,box_z,gamma,k_lr,vx,Nt)

      It will not make trajectory and output text files.


Version 0.1.7:
    - Now, gamma and k_lr is input parameters for the run_binding function.
      (not global constant)
      "run_binding(Nl,Nr,box_x,box_y,box_z,gamma,k_lr,vx,Nt,freq,o_name)"
      It will be helpful to adjust these values.

    - Now, make_acf function doesn't have struct output.
      just simply run the function. e.g.,
      b = make_acf(a.output_bit, cutoff), then b is the ACF array

    - make_xaxis function for making X-axis of ACF
      x = make_xaxis(cutoff, dt)
      then you can plot like this, plot(x, b)

Version 0.1.6:
    - Autocorrelation function added
      e.g., b = make_acf(a.output_bit, cutoff)
            you can access the acf in b.acf

    - New output
      a = run_binding(Nl, Nr, box_x, box_y, box_z, vx, Nt, freq, o_name)
      a.output:     raw data
      a.output_cut: cut the initial 1/10 trajectory of the raw data
      a.output_bit: if the count is larger than 2, the count is replaced by 1
                    (using output_cut)
      a.cumsum_bit: cumulative sum of output_bit

    - Simple codes for making initial positions
    - Pre-allocation for the update_pos_lig function
    - Normal distribution: changed to Random.randn
    - Some minor corrections

Version 0.1.3:
    - Some minor corrections

Version 0.1.2:
    - Now, output arrays (counts) will be saved to a text file.
    - Now, output arrays is 2D array to count binding events for each receptor

    - Added write_trj! function that can save trajectories.

    - Some minor corrections

Version 0.1.0:
    - Converting Pseudo-code to Julia

    - Note:
      Julia arrays are stored in column-major order,
      so pos_lig[3, Nl] array is about 100 times faster than
      pos_lig[Nl, 3] array for updating the force and position

    - To obtain a unit vector for force and position,
      abs() function was removed in dx, dy, and dz in update_force_output!

    - Some corrections for the reflecting boundary conditions of Y and Z
      e.g., if pos_lig[2,k] > D_box[2],
            then pos_lig[2,k] -= 2*(pos_lig[2,k]-D_box[2])
    - See the details in the update_pos_lig! function

    - So far, the global variables are just arbitrary values.

Pre-required module:
    - ProgressMeter (to see the progress)
    - DelimitedEfiles (to write output files)

Julia version:
    - 1.4.2 =#

##########################################################
# Intstuction

#= 
1. Import
    include("binding_v0.x.jl")
    using .LRBinding

2. Run
    a = run_binding(Nl, Nr, box_x, box_y, box_z, vx, Nt, freq, o_name)
    e.g., for 100 ligands, 2 receptors, [100,20,20] box, 0.01 velocity,
          10000 steps, writing a trajectory every 100 step,
          and "case_01.txt" and "case_01.xyz" output files
    a = run_binding(100, 2, 100, 20, 20, 0.01, 10000, 100, "case_01")

3. output
    One can access output array in a method.
    e.g., a.output

4. output.txt file is the transpose of output array =#

##########################################################
# Module LBBinding

module LRBinding

using ProgressMeter, DelimitedFiles, Random, StatsBase, LinearAlgebra
export run_binding, make_acf_dot, make_xaxis, analyze_results, make_acf3,
    analyze_results_acf3, calc_wait, trj2seg, seg2decimal, decimal2freq,
    sampleSn, calcLLF, trj2moreseg
# global variable (unmutable)

# dt = time step
const dt = 0.01

# The following attractive potential is removed in this code so there is no k_lr
# Non-bonding interactions between R and L, purely attractive

# U = ┬ k_lr * dis - a (dis <= d_rec/2), (k_lr > 0)
#     └ 0              (dis >  d_rec/2)

# F = -∇U = ┬ -k_lr    (dis <= d_rec/2)
#           └ 0        (dis >  d_rec/2)

# hence, -k_lr is constant force, and d_rec/2 is the cut-of

# Non-bonding interactions between L and L, WCA potential, purely repulsive

# U = ┬ 4*ep*((sigma/dis)^12 - (sigma/dis)^6) + ep        (dis <= 2^1/6*sigma)
#     └ 0                                                 (dis >  2^1/6*sigma)

# F = -∇U = ┬ 4*ep*(12*sigma^12/dis^13 - 6*sigma^6/dis^7) (dis <= 2^1/6*sigma)
#           └ 0                                           (dis >  2^1/6*sigma)

# hence, the cut-off distance is 2^(1/6)*sigma
# d_lig = 2^(1/6)*sigma * 2  # (2^(1/6) ~ 1.1246)
const epsilon, sigma = 0.01, 1.5

# mutable variable for output
mutable struct Results
    output::Vector{Int}
    pos_x::Vector{Float64}
    pos_y::Vector{Float64}
    pos_z::Vector{Float64}
end

mutable struct Results_ls
    freq::Vector{Int}
end

mutable struct Results_lsVec
    freq::Dict{Int,Vector{Int}}
end


mutable struct ACF
    acf::Vector{Float64}
    x::Vector{Float64}
    cdata::Vector{Float64}
    avg_btime::Float64
end

mutable struct Wtime
    wt_b_u::Vector{Int}
    wt_u_b::Vector{Int}
    wtime_output::Vector{Float64}
end
##########################################################
# Run

"""
This function is to execute the code.
It will produce an output as a method.

Nl: the number of ligand, Int
Nr: the number of receptor, Int
box_x, box_y, box_z: dimension of box, Real
gamma:friction coefficient, Float64
Temp: temperature, Float64
vx: receptor velocity, Float64
cutoff: cut-off radius for receptor, Float64
Eon: Energy of the binding state, Float64
Eoff: Energy of the unbinding state, FLoat64
Eb: Energy of the barrier, FLoat64
Nt: total time steps, Int
ls: length of segment for MLE
nth: sampling frequency
n: Sample size
freq: save trajectory every freq time step, Int
      if freq == 0, no trajectory.
o_name: file name for output(counts).txt and trajectory.xyz, String
"""
function run_binding(Nl::Int, Nr::Int, box_x::Real, box_y::Real, box_z::Real,
    gamma::Float64, alpha::Float64, Temp::Float64, vx::Float64,
    cutoff::Float64, Eon::Float64, Eoff::Float64, Eb::Float64,
    Nt::Int, freq::Int, o_name::String)
    output = zeros(Int, Nt) # counts
    D_box = Float64[box_x, box_y, box_z] # box dimension
    pos_x = D_box[1] .* rand(Nl)
    pos_y = D_box[2] .* rand(Nl)
    pos_z = D_box[3] .* rand(Nl)
    pos_rec = make_init_rec(D_box) # initial random position of receptor
    noise_x = Vector{Float64}(undef, Nl) # pre-allocation of noise X
    noise_y = Vector{Float64}(undef, Nl) # pre-allocation of noise Y
    noise_z = Vector{Float64}(undef, Nl) # pre-allocation of noise Z
    idx_lig = [0] # 0 for unbound receptor, or
    # index of a ligand that is bound to the receptor
    # equilibrium with 1 million steps
    println("1. Equilibrium step")
    p = Progress(1000000, 1) # for progress bar
    for i = 1:1000000
        update_idx_output!(pos_x, pos_y, pos_z, pos_rec, D_box, cutoff, Eon,
            Eoff, Eb, vx, gamma, alpha, idx_lig, Temp) # no output
        update_pos_lig!(pos_x, pos_y, pos_z, D_box, noise_x, noise_y, noise_z,
            gamma, Temp)
        update_pos_rec_idx!(pos_rec, vx, D_box, idx_lig, pos_x, pos_y, pos_z)
        next!(p) # for progress bar
    end
    # production run
    println("2. Production run")
    io = open(string(o_name, ".xyz"), "w") # trajectory file
    p = Progress(Nt, 1) # for progress bar
    for i = 1:Nt # iteration
        update_idx_output!(output, pos_x, pos_y, pos_z, pos_rec, D_box, cutoff,
            Eon, Eoff, Eb, vx, gamma, alpha, i, idx_lig, Temp)
        update_pos_lig!(pos_x, pos_y, pos_z, D_box, noise_x, noise_y, noise_z,
            gamma, Temp)
        update_pos_rec_idx!(pos_rec, vx, D_box, idx_lig, pos_x, pos_y, pos_z)
        if freq == 0 # no trajectory
        elseif rem(i, freq) == 0 # ever freq step
            write_trj!(pos_x, pos_y, pos_z, pos_rec, io)
        end
        next!(p) # for progress bar
    end
    close(io)
    if freq == 0 # remove the blank trajectory
        rm(string(o_name, ".xyz"))
    end
    # output(Nr, Nt) -> output(Nt, Nr), transpose for text file
    open(io -> writedlm(io, transpose(output)), string(o_name, ".txt"), "w")
    return Results(output, pos_x, pos_y, pos_z)
end

"""
A multiple dispatch of the run_binding function.
if one doesn't need to make output files of trajectory and output.txt,
just execute the run_binding function without the freq and o_name.
"""
function run_binding(Nl::Int, Nr::Int, box_x::Real, box_y::Real, box_z::Real,
    gamma::Float64, alpha::Float64, Temp::Float64, vx::Float64,
    cutoff::Float64, Eon::Float64, Eoff::Float64, Eb::Float64,
    Nt::Int)
    output = zeros(Int, Nt) # counts
    D_box = Float64[box_x, box_y, box_z] # box dimension
    # pos_lig = make_init_lig(Nl, D_box) # initial random position of ligand
    pos_x = D_box[1] .* rand(Nl)
    pos_y = D_box[2] .* rand(Nl)
    pos_z = D_box[3] .* rand(Nl)
    pos_rec = make_init_rec(D_box) # initial random position of receptor
    noise_x = Vector{Float64}(undef, Nl) # pre-allocation of noise X
    noise_y = Vector{Float64}(undef, Nl) # pre-allocation of noise Y
    noise_z = Vector{Float64}(undef, Nl) # pre-allocation of noise Z
    idx_lig = [0] # index of a ligand that is bound to the receptor, if the
    # receptor is in the unbinding state, then 0
    # equilibrium with 1 million steps
    println("1. Equilibrium step")
    p = Progress(1000000, 1) # for progress bar
    for i = 1:1000000
        update_idx_output!(pos_x, pos_y, pos_z, pos_rec, D_box, cutoff, Eon,
            Eoff, Eb, vx, gamma, alpha, idx_lig, Temp) # no output
        update_pos_lig!(pos_x, pos_y, pos_z, D_box, noise_x, noise_y, noise_z,
            gamma, Temp)
        update_pos_rec_idx!(pos_rec, vx, D_box, idx_lig, pos_x, pos_y, pos_z)
        next!(p) # for progress bar
    end
    # production run
    println("2. Production run")
    p = Progress(Nt, 1) # for progress bar
    for i = 1:Nt # iteration
        update_idx_output!(output, pos_x, pos_y, pos_z, pos_rec, D_box, cutoff,
            Eon, Eoff, Eb, vx, gamma, alpha, i, idx_lig, Temp)
        update_pos_lig!(pos_x, pos_y, pos_z, D_box, noise_x, noise_y, noise_z,
            gamma, Temp)
        update_pos_rec_idx!(pos_rec, vx, D_box, idx_lig, pos_x, pos_y, pos_z)
        next!(p) # for progress bar
    end
    return Results(output, pos_x, pos_y, pos_z)
end

"""
A multiple dispatch of the run_binding function.
It won't give us output but probability for ls and nth.
ls and nth are an Int type.
for example, ls = 2
freq = [P00, P01, P10, P11])
"""
function run_binding(Nl::Int, Nr::Int, box_x::Real, box_y::Real, box_z::Real,
    gamma::Float64, alpha::Float64, Temp::Float64, vx::Float64,
    cutoff::Float64, Eon::Float64, Eoff::Float64, Eb::Float64,
    Nt::Int, ls::Int, nth::Int)
    #output = zeros(Int, Nt) # counts
    D_box = Float64[box_x, box_y, box_z] # box dimension
    pos_x = D_box[1] .* rand(Nl)
    pos_y = D_box[2] .* rand(Nl)
    pos_z = D_box[3] .* rand(Nl)
    pos_rec = reshape(D_box / 2, 3, 1)
    noise_x = Vector{Float64}(undef, Nl) # pre-allocation of noise X
    noise_y = Vector{Float64}(undef, Nl) # pre-allocation of noise Y
    noise_z = Vector{Float64}(undef, Nl) # pre-allocation of noise Z
    idx_lig = [0] # index of a ligand that is bound to the receptor, if the
    # receptor is in the unbinding state, then 0
    seg = zeros(Int, nth, ls) # pre-allocation of seg
    freq = zeros(Int, 2^ls) # output
    # equilibrium with 1 million steps
    println("1. Equilibrium step")
    p = Progress(1000000, 1) # for progress bar
    for i = 1:1000000
        update_idx_output!(pos_x, pos_y, pos_z, pos_rec, D_box, cutoff, Eon,
            Eoff, Eb, vx, gamma, alpha, idx_lig, Temp) # no output
        update_pos_lig!(pos_x, pos_y, pos_z, D_box, noise_x, noise_y, noise_z,
            gamma, Temp)
        update_pos_rec_idx!(pos_rec, vx, D_box, idx_lig, pos_x, pos_y, pos_z)
        next!(p) # for progress bar
    end
    # production run
    println("2. Production run")
    p = Progress(Nt, 1) # for progress bar
    for i = 1:Nt # iteration
        output = update_idx_output(pos_x, pos_y, pos_z, pos_rec, D_box, cutoff,
            Eon, Eoff, Eb, vx, gamma, alpha, idx_lig, Temp)
        countFreq!(seg, freq, i, output, ls, nth)
        update_pos_lig!(pos_x, pos_y, pos_z, D_box, noise_x, noise_y, noise_z,
            gamma, Temp)
        update_pos_rec_idx!(pos_rec, vx, D_box, idx_lig, pos_x, pos_y, pos_z)
        next!(p) # for progress bar
    end
    return Results_ls(freq)
end

"""
A multiple dispatch of the run_binding function.
It won't give us output but probability for ls and nth.
lsVec a vector that contains ls, and nth are an Int type.
for example, ls = [2,3]
freq = Dict[2 => [P00, P01, P10, P11]
            3 => [P000, P001, P010, P011, P100, P101, P110, P111]]
"""
function run_binding(Nl::Int, Nr::Int, box_x::Real, box_y::Real, box_z::Real,
    gamma::Float64, alpha::Float64, Temp::Float64, vx::Float64,
    cutoff::Float64, Eon::Float64, Eoff::Float64, Eb::Float64,
    Nt::Int, lsVec::Vector{Int}, nth::Int)
    #output = zeros(Int, Nt) # counts
    D_box = Float64[box_x, box_y, box_z] # box dimension
    pos_x = D_box[1] .* rand(Nl)
    pos_y = D_box[2] .* rand(Nl)
    pos_z = D_box[3] .* rand(Nl)
    pos_rec = reshape(D_box / 2, 3, 1)
    noise_x = Vector{Float64}(undef, Nl) # pre-allocation of noise X
    noise_y = Vector{Float64}(undef, Nl) # pre-allocation of noise Y
    noise_z = Vector{Float64}(undef, Nl) # pre-allocation of noise Z
    idx_lig = [0] # index of a ligand that is bound to the receptor, if the
    # receptor is in the unbinding state, then 0
    seg = Dict(ls => zeros(Int, nth, ls) for ls in lsVec)
    freq = Dict(ls => zeros(Int, 2^ls) for ls in lsVec) # output
    # equilibrium with 1 million steps
    println("1. Equilibrium step")
    p = Progress(1000000, 1) # for progress bar
    for i = 1:1000000
        update_idx_output!(pos_x, pos_y, pos_z, pos_rec, D_box, cutoff, Eon,
            Eoff, Eb, vx, gamma, alpha, idx_lig, Temp) # no output
        update_pos_lig!(pos_x, pos_y, pos_z, D_box, noise_x, noise_y, noise_z,
            gamma, Temp)
        update_pos_rec_idx!(pos_rec, vx, D_box, idx_lig, pos_x, pos_y, pos_z)
        next!(p) # for progress bar
    end
    # production run
    println("2. Production run")
    p = Progress(Nt, 1) # for progress bar
    for i = 1:Nt # iteration
        output = update_idx_output(pos_x, pos_y, pos_z, pos_rec, D_box, cutoff,
            Eon, Eoff, Eb, vx, gamma, alpha, idx_lig, Temp)
        countFreq!(seg, freq, i, output, nth)
        update_pos_lig!(pos_x, pos_y, pos_z, D_box, noise_x, noise_y, noise_z,
            gamma, Temp)
        update_pos_rec_idx!(pos_rec, vx, D_box, idx_lig, pos_x, pos_y, pos_z)
        next!(p) # for progress bar
    end
    return Results_lsVec(freq)
end


##########################################################
# Core functions
""" make an initial position of ligands """
function make_init_lig(Nl::Int, D_box::Vector{Float64})
    return D_box .* rand(3, Nl)
end

""" make an initial position of receptors """
function make_init_rec(D_box::Vector{Float64})
    return reshape(D_box / 2, 3, 1)
end

function apply_boundary!(pos_x::Vector{Float64}, pos_y::Vector{Float64},
    pos_z::Vector{Float64}, D_box::Vector{Float64})
    # Apply the reflecting boundary and PBC
    @inbounds for k in eachindex(pos_x)
        # PBC x
        if pos_x[k] > D_box[1]
            pos_x[k] -= D_box[1]
        elseif pos_x[k] <= 0
            pos_x[k] += D_box[1]
        end
        # reflecting Y
        if pos_y[k] > D_box[2]
            pos_y[k] -= 2 * (pos_y[k] - D_box[2])
        elseif pos_y[k] <= 0
            pos_y[k] = -pos_y[k]
        end
        # reflecting Z
        if pos_z[k] > D_box[3]
            pos_z[k] -= 2 * (pos_z[k] - D_box[3])
        elseif pos_z[k] <= 0
            pos_z[k] = -pos_z[k]
        end
    end
    return nothing
end

""" update idx_ligand and output """
function update_idx_output!(output::Vector{Int}, pos_x::Vector{Float64},
    pos_y::Vector{Float64}, pos_z::Vector{Float64}, pos_rec::Array{Float64,2},
    D_box::Vector{Float64}, cutoff::Float64, Eon::Float64, Eoff::Float64,
    Eb::Float64, vx::Float64, gamma::Float64, alpha::Float64, i::Int,
    idx_lig::Vector{Int}, Temp::Float64)
    within_lig = Int.([])
    if idx_lig[1] == 0
        # finding ligands within cutoff radius
        @inbounds for k in eachindex(pos_x) # 1:Nl
            dx = pos_x[k] - pos_rec[1, 1]
            dx -= D_box[1] * round(dx / D_box[1])  # for PBC
            dy = pos_y[k] - pos_rec[2, 1]
            dz = pos_z[k] - pos_rec[3, 1]
            if abs(dx) <= cutoff && abs(dy) <= cutoff && abs(dz) <= cutoff
                dis = sqrt(dx^2 + dy^2 + dz^2)
                if dis <= cutoff
                    push!(within_lig, k)
                end
            end
        end
        if length(within_lig) == 0
            output[i] = 0
        else
            P_on = length(within_lig) * dt * exp(-(Eb - Eoff) / Temp)
            if rand() <= P_on
                idx_lig[1] = rand(within_lig)
                output[i] = 1
                pos_x[idx_lig[1]] = pos_rec[1, 1]
                pos_y[idx_lig[1]] = pos_rec[2, 1]
                pos_z[idx_lig[1]] = pos_rec[3, 1]
            else
                output[i] = 0
            end
        end
    else
        P_off = dt * exp(-(Eb - Eon - (gamma * alpha * vx)) / Temp)
        if rand() <= P_off
            idx_lig[1] = 0
            output[i] = 0
        else
            output[i] = 1
        end
    end
    return nothing
end

""" update idx_ligand without output for equilibrium step """
function update_idx_output!(pos_x::Vector{Float64}, pos_y::Vector{Float64},
    pos_z::Vector{Float64}, pos_rec::Array{Float64,2}, D_box::Vector{Float64},
    cutoff::Float64, Eon::Float64, Eoff::Float64, Eb::Float64, vx::Float64,
    gamma::Float64, alpha::Float64, idx_lig::Vector{Int}, Temp::Float64)
    within_lig = Int.([])
    if idx_lig[1] == 0
        # finding ligands within cutoff radius
        @inbounds for k in eachindex(pos_x) # 1:Nl
            dx = pos_x[k] - pos_rec[1, 1]
            dx -= D_box[1] * round(dx / D_box[1])  # for PBC
            dy = pos_y[k] - pos_rec[2, 1]
            dz = pos_z[k] - pos_rec[3, 1]
            if abs(dx) <= cutoff && abs(dy) <= cutoff && abs(dz) <= cutoff
                dis = sqrt(dx^2 + dy^2 + dz^2)
                if dis <= cutoff
                    push!(within_lig, k)
                end
            end
        end
        if length(within_lig) != 0
            P_on = length(within_lig) * dt * exp(-(Eb - Eoff) / Temp)
            if rand() <= P_on
                idx_lig[1] = rand(within_lig)
                # output[i] = 1
                pos_x[idx_lig[1]] = pos_rec[1, 1]
                pos_y[idx_lig[1]] = pos_rec[2, 1]
                pos_z[idx_lig[1]] = pos_rec[3, 1]
            end
        end
    else
        P_off = dt * exp(-(Eb - Eon - (gamma * alpha * vx)) / Temp)
        if rand() <= P_off
            idx_lig[1] = 0
            # output[i] = 0
            # else
            #    output[i] = 1
        end
    end
    return nothing
end

""" update idx_ligand and output for freq calculation"""
function update_idx_output(pos_x::Vector{Float64}, pos_y::Vector{Float64},
    pos_z::Vector{Float64}, pos_rec::Array{Float64,2}, D_box::Vector{Float64},
    cutoff::Float64, Eon::Float64, Eoff::Float64, Eb::Float64, vx::Float64,
    gamma::Float64, alpha::Float64, idx_lig::Vector{Int}, Temp::Float64)
    within_lig = Int.([])
    output = 0
    if idx_lig[1] == 0
        # finding ligands within cutoff radius
        @inbounds for k in eachindex(pos_x) # 1:Nl
            dx = pos_x[k] - pos_rec[1, 1]
            dx -= D_box[1] * round(dx / D_box[1])  # for PBC
            dy = pos_y[k] - pos_rec[2, 1]
            dz = pos_z[k] - pos_rec[3, 1]
            if abs(dx) <= cutoff && abs(dy) <= cutoff && abs(dz) <= cutoff
                dis = sqrt(dx^2 + dy^2 + dz^2)
                if dis <= cutoff
                    push!(within_lig, k)
                end
            end
        end
        if length(within_lig) != 0
            P_on = length(within_lig) * dt * exp(-(Eb - Eoff) / Temp)
            if rand() <= P_on
                idx_lig[1] = rand(within_lig)
                output = 1
                pos_x[idx_lig[1]] = pos_rec[1, 1]
                pos_y[idx_lig[1]] = pos_rec[2, 1]
                pos_z[idx_lig[1]] = pos_rec[3, 1]
            end
        end
    else
        P_off = dt * exp(-(Eb - Eon - (gamma * alpha * vx)) / Temp)
        if rand() <= P_off
            idx_lig[1] = 0
        else
            output = 1
        end
    end
    return output
end

""" get index """
function getIdx(i::Int, ls::Int, nth::Int)
    idx = rem(i, ls * nth)
    if idx == 0
        idx = ls * nth
    end
    return idx
end

function getCoord(idx::Int, ls::Int, nth::Int)
    y, x = divrem(idx, nth)
    y += 1
    if x == 0
        y -= 1
        x = nth
    end
    return x, y
end

function getDecimal(seg::Array{Int,2}, x::Int, y::Int, ls::Int)
    decimal = 0
    for i = 1:ls
        y += 1
        if y > ls
            y = 1
        end
        decimal += seg[x, y] * 2^(ls - i)
    end
    decimal += 1
    return decimal
end

function countFreq!(seg::Array{Int,2}, freq::Vector{Int}, i::Int, output::Int,
    ls::Int, nth::Int)
    idx = getIdx(i, ls, nth)
    seg[idx] = output
    x, y = getCoord(idx, ls, nth)
    decimal = getDecimal(seg, x, y, ls)
    if i > (ls - 1) * nth
        freq[decimal] += 1
    end
    return nothing
end

function countFreq!(seg::Dict{Int,Array{Int,2}}, freq::Dict{Int,Vector{Int}},
    i::Int, output::Int, nth::Int)
    for key in keys(seg)
        idx = getIdx(i, key, nth)
        seg[key][idx] = output
        x, y = getCoord(idx, key, nth)
        decimal = getDecimal(seg[key], x, y, key)
        if i > (key - 1) * nth
            freq[key][decimal] += 1
        end
    end
    return nothing
end

""" update ligand positions """
function update_pos_lig!(pos_x::Vector{Float64}, pos_y::Vector{Float64},
    pos_z::Vector{Float64}, D_box::Vector{Float64}, noise_x::Vector{Float64},
    noise_y::Vector{Float64}, noise_z::Vector{Float64}, gamma::Float64,
    Temp::Float64)
    randn!(noise_x)
    randn!(noise_y)
    randn!(noise_z)
    noise_k = sqrt(2 * Temp * dt / gamma)
    # update positions using the Brownian dynamics
    @inbounds for k in eachindex(pos_x)
        pos_x[k] += noise_k * noise_x[k]
        pos_y[k] += noise_k * noise_y[k]
        pos_z[k] += noise_k * noise_z[k]
    end
    # Apply the reflecting boundary and PBC
    apply_boundary!(pos_x, pos_y, pos_z, D_box)
    return nothing
end

""" update receptor and bound ligand(if present) positions """
function update_pos_rec_idx!(pos_rec::Array{Float64,2}, vx::Float64,
    D_box::Vector{Float64}, idx_lig::Vector{Int}, pos_x::Vector{Float64},
    pos_y::Vector{Float64}, pos_z::Vector{Float64})
    # only along X direction
    pos_rec[1, 1] += dt * vx
    # PBC x
    if pos_rec[1, 1] > D_box[1]
        pos_rec[1, 1] -= D_box[1]
    elseif pos_rec[1, 1] <= 0
        pos_rec[1, 1] += D_box[1]
    end
    if idx_lig[1] != 0 # update the position of the ligand bound to the receptor
        pos_x[idx_lig[1]] = pos_rec[1, 1]
        pos_y[idx_lig[1]] = pos_rec[2, 1]
        pos_z[idx_lig[1]] = pos_rec[3, 1]
    end
    return nothing
end

""" write trajectory file to XYZ file format """
function write_trj!(pos_x::Vector{Float64}, pos_y::Vector{Float64},
    pos_z::Vector{Float64}, pos_rec::Array{Float64,2}, io::IOStream)
    n_particle = length(pos_x) + size(pos_rec, 2)
    write(io, string(n_particle), "\n\n") # number of particle + blank line
    for j in eachindex(pos_x)
        write(io, "C\t") # atom name C for ligand
        writedlm(io, round.([pos_x[j] pos_y[j] pos_z[j]], digits=3))
    end
    for j in axes(pos_rec, 2)
        write(io, "N\t") # atom name N for receptor
        writedlm(io, round.([pos_rec[1, j] pos_rec[2, j] pos_rec[3, j]], digits=3))
    end
    return nothing
end

###########################################################
# MLE functions

""" Converting trajectory to segment """
function trj2seg(output::Vector{Int}, ls::Int, nth::Int)
    # ls is the number of bits in the segment, 
    # nth*delT is the time difference between two adjacent bits
    Nt = length(output) #number \of total time steps
    ns = trunc(Int, Nt / ls / nth) #number of segments
    seg = zeros(Int, (ls, ns))
    for i = 1:ns
        for j = 1:ls
            seg[j, i] = output[(i-1)*nth*ls+(j-1)*nth+1]
        end
    end
    return seg
end

# this function is the same as trj2seg. No performance improvement, but simple.
""" Converting trajectory to segment """
function trj2seg2(output::Vector{Int}, ls::Int, nth::Int)
    seg = reshape(output[begin:nth:end], ls, :)
    return seg
end

function old_trj2moreseg(output::Vector{Int}, ls::Int, nth::Int)
    # ls is the number of bits in the segment, 
    # nth*delT is the time difference between two adjacent bits
    Nt = length(output) #number \of total time steps
    ns = Nt - nth #number of segments
    seg = zeros(Int, (ls, ns))
    for i = 1:ns
        seg[1, i] = output[i]
        seg[2, i] = output[i+nth]
    end
    return seg
end

function trj2moreseg(output::Vector{Int}, ls::Int, nth::Int)
    # ls is the number of bits in the segment, 
    # nth*delT is the time difference between two adjacent bits
    Nt = length(output) #number \of total time steps
    ns = Nt - nth * (ls - 1) #number of segments
    seg = zeros(Int, (ls, ns))
    for i = 1:ls
        seg[i, :] = output[(1+nth*(i-1)):(nth*(i-1)+ns)]
    end
    return seg
end

""" Converting segment array to decimal vector """
function seg2decimal(seg::Array{Int,2})
    # binary is a column vector of length ls
    binary = [2^(i) for i = 0:(size(seg, 1)-1)]
    decimal = (transpose(reverse(binary))*seg)[:] .+ 1
    return decimal
end

""" Collecting frequency for each decimal in the sample space """
function decimal2freq(decimal::Array{Int}, ls::Int)
    freq = zeros(Int, 2^ls)
    for i = eachindex(decimal)
        a = decimal[i]
        freq[a] += 1
    end
    return freq
end

""" sampling Sn from frequency """
function freq2Sn(freq::Vector{Int}, nsize::Int)
    cumsumRate = cumsum(freq) / sum(freq)
    Sn = Vector{Int}(undef, nsize)
    for i = 1:nsize
        Sn[i] = searchsortedfirst(cumsumRate, rand())
    end
    return Sn
end

""" old one, Calculating log likelihood function """
function old_calcLLF(Sn::Vector{Int}, freq::Vector{Int})
    llf = 0
    tot = sum(freq)
    for i = eachindex(Sn)
        a = Sn[i]
        llf += log(freq[a] / tot)
    end
    return llf
end

""" Calculating log likelihood function """
function calcLLF(freq_star::Vector{Int}, freq::Vector{Int})
    if sum(freq_star) != sum(freq)
        throw(DomainError("samples are not matching in calcLLF"))
    end
    llf = Float64(0)
    tot = sum(freq)
    for i = eachindex(freq_star)
        llf += freq_star[i] * log(freq[i])
    end
    llf = llf / tot - log(tot)
    return llf
end

function calcLLF_inside(freq_matrix::Array{Float64,2})
    hx = 2.0
    hy = 0.02
    hz = 0.002
    llf = Vector{Float64}(undef, 10)
    freq_star = freq_matrix[7, :]
    tot = sum(freq_star)
    for i = 1:size(freq_matrix, 1)
        llf_init = 0
        for j = 1:size(freq_matrix, 2)
            llf_init += freq_star[j] * log(freq_matrix[i, j])
        end
        llf_init = llf_init / tot - log(tot)
        llf[i] = llf_init
    end
    # hessian
    hessian = zeros(Float64, 3, 3)
    firstD = zeros(Float64, 3)
    firstD[1] = (llf[10] - llf[3]) / (2 * hx)
    firstD[2] = (llf[9] - llf[5]) / (2 * hy)
    firstD[3] = (llf[8] - llf[6]) / (2 * hz)
    hessian[1, 1] = (llf[3] + llf[10] - 2 * llf[7]) / (hx^2)
    hessian[2, 2] = (llf[5] + llf[9] - 2 * llf[7]) / (hy^2)
    hessian[3, 3] = (llf[6] + llf[8] - 2 * llf[7]) / (hz^2)

    hessian[1, 2] = (llf[1] - llf[5] - llf[3] + llf[7]) / (hx * hy)
    hessian[1, 3] = (llf[2] - llf[6] - llf[3] + llf[7]) / (hx * hz)
    hessian[2, 3] = (llf[4] - llf[6] - llf[5] + llf[7]) / (hz * hy)

    hessian[2, 1] = hessian[1, 2]
    hessian[3, 1] = hessian[1, 3]
    hessian[3, 2] = hessian[2, 3]
    return hessian, firstD
end



##########################################################
# Analysis

""" Analysis of ACF and average binding time"""
function analyze_results(output::Vector{Int}, lags::Int, n::Int)
    if lags * dt < n
        throw(DomainError("n exceeds the lags * dt."))
    end
    acf = make_acf_dot(output, lags) # ACF
    xdata = make_xaxis(lags, dt) # x-axis
    cdata = calc_cn(acf, dt, n) # C(n) vector
    avg_btime = sum(output) / length(output) # average binding time
    return ACF(acf, xdata, cdata, avg_btime)
end

""" Analysis of ACF3 and average binding time"""
function analyze_results_acf3(output::Vector{Int}, lags::Int, n::Int)
    if lags * dt < n
        throw(DomainError("n exceeds the lags * dt."))
    end
    acf = make_acf3(output, lags) # ACF
    xdata = make_xaxis(lags, dt) # x-axis
    cdata = calc_cn(acf, dt, n) # C(n) vector
    avg_btime = sum(output) / length(output) # average binding time
    return ACF(acf, xdata, cdata, avg_btime)
end

"""
Autocorrelation function using dot product
Input parameter (output) is a vector
"""
function make_acf_dot(output::Vector{Int}, lags::Int)
    nstep = length(output)
    acf_temp = Vector{Float64}(undef, lags + 1)
    for j = 0:lags
        acf_temp[j+1] = dot(view(output, 1:(lastindex(output)-j)),
            view(output, (1+j):lastindex(output)))
    end
    count_num = collect(nstep:-1:(nstep-lags))
    acf = acf_temp ./ count_num # mean
    acf_norm = acf / acf[1] # normalized
    return acf_norm
end

"""
triple autocorrelation function C3, especially f(x)f(x+τ)f(x+2τ)
Input parameter (output) is a vector
"""
function make_acf3(output::Vector{Int}, lags::Int)
    nstep = length(output)
    acf_temp = Vector{Float64}(undef, lags + 1)
    for j = 0:lags
        h1 = view(output, 1:(lastindex(output)-2*j))
        h2 = view(output, (1+j):lastindex(output)-j)
        h3 = view(output, (1+2*j):lastindex(output))
        acf_temp[j+1] = elem_multi_3v(h1, h2, h3)
    end
    count_num = collect(nstep:-2:(nstep-2*lags))
    acf = acf_temp ./ count_num # mean
    acf_norm = acf / acf[1] # normalized
    return acf_norm
end

""" sub function element-wise multiplication for three vectors """
function elem_multi_3v(a::AbstractVector{Int}, b::AbstractVector{Int},
    c::AbstractVector{Int})
    if !(length(a) == length(b) == length(c))
        throw(DomainError("the length of three vectors is not the same"))
    end
    return sum(Ia * Ib * Ic for (Ia, Ib, Ic) in zip(a, b, c))
end

""" X-axis for ACF """
function make_xaxis(lags::Int, dt::Float64)
    return collect(0:dt:lags*dt)
end

""" Calculate C(n) """
function calc_cn(acf::Vector{Float64}, dt::Float64, n::Int)
    cdata = Vector{Float64}(undef, n) # C(1), C(2), ... ,C(n)
    for i = 1:n
        cdata[i] = acf[1+convert(Int, round(i / dt))]
    end
    return cdata
end

""" wating time (input parameter is a vector) """
function calc_wait(output::Vector{Int})
    wt_b_u = Vector{Int}() # waiting time from binding to unbinding
    wt_u_b = Vector{Int}() # waiting time from unbinding to binding
    wtime_output = Vector{Float64}(undef, 6)
    wt_whole = rle(output) # run-length encoding of a vector as a tuple
    if wt_whole[1][begin] == 1 # if trajectory starts from binding
        append!(wt_b_u, wt_whole[2][(begin+2):2:end] .- 1) # remove first
        append!(wt_u_b, wt_whole[2][(begin+1):2:end] .- 1)
    else # if trajectory starts from unbinding
        append!(wt_b_u, wt_whole[2][(begin+1):2:end] .- 1)
        append!(wt_u_b, wt_whole[2][(begin+2):2:end] .- 1) # remove first
    end
    if wt_whole[1][end] == 1 # remove the last one
        pop!(wt_b_u)
    else
        pop!(wt_u_b)
    end
    # wtime_u = fit(Histogram, wt_b_u, minimum(wt_b_u):maximum(wt_b_u) + 1)
    # wtime_b = fit(Histogram, wt_u_b, minimum(wt_u_b):maximum(wt_u_b) + 1)
    m1_b_u = mean(wt_b_u) # first moment (b -> u)
    m2_b_u = mean(wt_b_u .^ 2) # second moment (b -> u)
    m3_b_u = mean(wt_b_u .^ 3) # third moment (b -> u)
    m1_u_b = mean(wt_u_b) # first moment (u -> b)
    m2_u_b = mean(wt_u_b .^ 2) # second moment (u -> b)
    m3_u_b = mean(wt_u_b .^ 3) # third moment (u -> b)
    wtime_output = [m1_b_u; m2_b_u; m3_b_u; m1_u_b; m2_u_b; m3_u_b]
    return Wtime(wt_b_u, wt_u_b, wtime_output)
end

##############################################################

function calc_energy(pos_x::Vector{Float64}, pos_y::Vector{Float64},
    pos_z::Vector{Float64}, D_box::Vector{Float64}, kloop::Vector{Int},
    jloop::Vector{Int})
    # initialize
    cutd = 2^(1 / 6) * sigma
    energy_x = zeros(Float64, length(pos_x))
    energy_y = zeros(Float64, length(pos_x))
    energy_z = zeros(Float64, length(pos_x))
    @inbounds for l in eachindex(kloop)
        dx = pos_x[kloop[l]] - pos_x[jloop[l]]
        dx -= D_box[1] * round(dx / D_box[1])
        dy = pos_y[kloop[l]] - pos_y[jloop[l]]
        dz = pos_z[kloop[l]] - pos_z[jloop[l]]
        # d_lig/2 = 2^(1/6)*sigma
        if abs(dx) <= cutd && abs(dy) <= cutd && abs(dz) <= cutd
            dis = sqrt(dx^2 + dy^2 + dz^2)
            if dis <= cutd
                abs_energy = 4 * epsilon * (sigma^12 / dis^12 -
                                            sigma^6 / dis^6) + epsilon
                ex = abs_energy * dx / dis
                ey = abs_energy * dy / dis
                ez = abs_energy * dz / dis
                energy_x[kloop[l]] += ex
                energy_y[kloop[l]] += ey
                energy_z[kloop[l]] += ez
                energy_x[jloop[l]] -= ex
                energy_y[jloop[l]] -= ey
                energy_z[jloop[l]] -= ez
            end
        end
    end
    return energy_x, energy_y, energy_z
end

""" update force """
function update_force!(force_x::Vector{Float64}, force_y::Vector{Float64},
    force_z::Vector{Float64}, pos_x::Vector{Float64}, pos_y::Vector{Float64},
    pos_z::Vector{Float64}, D_box::Vector{Float64}, kloop::Vector{Int},
    jloop::Vector{Int})
    # initialize
    cutd = 2^(1 / 6) * sigma
    fill!(force_x, 0.0)
    fill!(force_y, 0.0)
    fill!(force_z, 0.0)
    @inbounds for l in eachindex(kloop) # 1:Nl
        dx = pos_x[kloop[l]] - pos_x[jloop[l]]
        dx -= D_box[1] * round(dx / D_box[1])
        dy = pos_y[kloop[l]] - pos_y[jloop[l]]
        dz = pos_z[kloop[l]] - pos_z[jloop[l]]
        # d_lig/2 = 2^(1/6)*sigma
        if abs(dx) <= cutd && abs(dy) <= cutd && abs(dz) <= cutd
            dis = sqrt(dx^2 + dy^2 + dz^2)
            if dis <= cutd
                abs_force = 4 * epsilon * (12 * sigma^12 / dis^13 -
                                           6 * sigma^6 / dis^7)
                fx = abs_force * dx / dis
                fy = abs_force * dy / dis
                fz = abs_force * dz / dis
                force_x[kloop[l]] += fx
                force_y[kloop[l]] += fy
                force_z[kloop[l]] += fz
                force_x[jloop[l]] -= fx
                force_y[jloop[l]] -= fy
                force_z[jloop[l]] -= fz
            end
        end
    end
    return nothing
end

function new_position!(new_position::AbstractVector{Float64},
    pos::AbstractVector{Float64}, force::AbstractVector{Float64},
    max_force::Float64, h::Float64)
    @inbounds for i in eachindex(pos)
        new_position[i] = pos[i] + force[i] * h / max_force
    end
    return nothing
end

function idx_energy_accept!(idx::BitVector, energy::Vector{Float64},
    new_energy::Vector{Float64})
    @inbounds for i in eachindex(energy)
        idx[i] = new_energy[i] < energy[i]
    end
    return nothing
end

""" Energy minimization """
function steepest_descent!(pos_x::Vector{Float64}, pos_y::Vector{Float64},
    pos_z::Vector{Float64}, force_x::Vector{Float64}, force_y::Vector{Float64},
    force_z::Vector{Float64}, D_box::Vector{Float64}, kloop::Vector{Int},
    jloop::Vector{Int}, nstep::Int)
    h = fill(0.1, length(pos_x)) # initial h value, 0.1 is enough for our system
    for i = 1:nstep
        update_force!(force_x, force_y, force_z, pos_x, pos_y, pos_z, D_box,
            kloop, jloop)
        maxF = maximum(sqrt.(force_x .^ 2 + force_y .^ 2 + force_z .^ 2))
        if maxF == 0.0
            println("Energy minimized at $i steps")
            return nothing
        end
        energy_x, energy_y, energy_z = calc_energy(pos_x, pos_y, pos_z, D_box,
            kloop, jloop)
        new_pos_x = pos_x .+ force_x ./ maxF .* h
        new_pos_y = pos_y .+ force_y ./ maxF .* h
        new_pos_z = pos_z .+ force_z ./ maxF .* h
        apply_boundary!(new_pos_x, new_pos_y, new_pos_z, D_box)
        new_energy_x, new_energy_y, new_energy_z = calc_energy(new_pos_x,
            new_pos_y, new_pos_z, D_box, kloop, jloop) # new energy calc
        idx_x = new_energy_x .< energy_x
        idx_y = new_energy_y .< energy_y
        idx_z = new_energy_z .< energy_z
        pos_x[idx_x] = new_pos_x[idx_x] # accept
        pos_y[idx_y] = new_pos_y[idx_y] # accept
        pos_z[idx_z] = new_pos_z[idx_z] # accept
    end
    return nothing
end

function calc_closest_lig(pos_x::Vector{Float64}, pos_y::Vector{Float64},
    pos_z::Vector{Float64}, pos_rec::Array{Float64,2}, within_lig::Vector{Int},
    D_box::Vector{Float64})
    closest = Int(0)
    temp = Float64(0)
    for i = eachindex(within_lig)
        k = within_lig[i]
        dx = pos_x[k] - pos_rec[1, 1]
        dx -= D_box[1] * round(dx / D_box[1])  # for PBC
        dy = pos_y[k] - pos_rec[2, 1]
        dz = pos_z[k] - pos_rec[3, 1]
        dis = sqrt(dx^2 + dy^2 + dz^2)
        if closest == 0
            closest = k
            temp = dis
        elseif dis < temp
            closest = k
            temp = dis
        end
    end
    return closest
end

function kj_loop(Nl::Int)
    kloop = Int[]
    jloop = Int[]
    for k in 1:Nl
        for j in 1:Nl
            if j < k
                push!(kloop, k)
                push!(jloop, j)
            end
        end
    end
    return kloop, jloop
end


##########################################################
# Post analysis (combined)

##########################################################
# removed

#= 

mutable struct ACF
    converged::Vector{Bool}
    param::Array{Float64,2}
    fitting::Array{Float64,2}
    lifetime::Vector{Float64}
    x::Vector{Float64}
    plot::Array{Float64,2}
end

""" Autocorrelation function """
function make_acf(output_bit::Array{Int,2}, cutoff::Int)
    output_acf = Array{Float64,2}(undef, (cutoff + 1, size(output_bit, 1)))
    for i = axes(output_bit, 1)
        output_acf[:,i] = autocor(output_bit[i,:], 0:cutoff)
    end
    return output_acf
end

""" Autocorrelation function using for-loops """
function make_acf_forloop(output_bit::Array{Int,2}, lags::Int)
    nr, nstep = size(output_bit) # number of receptors and time steps
    acf = Array{Float64,2}(undef, (lags + 1, nr))
    temp = Vector{Int}(undef, size(output_bit, 2))
    for i = axes(output_bit, 1) # for each receptor
        copyto!(temp, view(output_bit, i, :)) # copy ith receptor to temp vector
        for j = 0:lags # from 0 (self) to the maximum lag
            sum = 0
            for k = 1:(nstep - j)
                sum += temp[k] * temp[k + j]
            end
            acf[j + 1,i] = sum / (nstep - j) # mean
        end
    end
    output_acf = acf ./ transpose(acf[1,:]) # normalized
    return output_acf
end

""" 3-term exponential fitting """
function threeExpFit(output_acf::Array{Float64,2}, xdata::Vector{Float64},
                     p0::Vector{Float64})
    model(t, p) = p[1] * exp.(-t / p[2]) + p[3] * exp.(-t / p[4]) +
                  (1 - p[1] - p[3]) * exp.(-t / p[5]) # 3-term exponential
    lb = [0.0, 0.0, 0.0, 0.0, 0.0] # lower bounds
    output_converged = Vector{Bool}(undef, size(output_acf, 2))
    output_param = Array{Float64,2}(undef, 5, size(output_acf, 2))
    output_fitting = Array{Float64,2}(undef, size(output_acf))
    output_lifetime = Vector{Float64}(undef, size(output_acf, 2))
    for i = axes(output_acf, 2)
        fit = curve_fit(model, xdata, output_acf[:,i], p0, lower=lb)
        output_converged[i] = fit.converged
        output_param[:,i] = fit.param
        output_fitting[:,i] = model(xdata, fit.param)
        output_lifetime[i] = fit.param[1] * fit.param[2] +
                             fit.param[3] * fit.param[4] +
                             (1 - fit.param[1] - fit.param[3]) * fit.param[5]
    end
    return output_converged, output_param, output_fitting, output_lifetime
end

""" Autocorrelation function using for dot product """
function make_acf_dot(output_bit::Array{Int,2}, lags::Int)
    nr, nstep = size(output_bit)
    acf = Array{Float64,2}(undef, (lags + 1, nr))
    temp = Vector{Int}(undef, nstep)
    for i = axes(output_bit, 1)
        temp[:] = view(output_bit, i, :) # copy
        for j = 0:lags
            acf[j + 1, i] = dot(view(temp, 1:(lastindex(temp) - j)),
                view(temp, (1 + j):lastindex(temp))) / (nstep - j)
        end
    end
    output_acf = acf ./ transpose(acf[1,:]) # normalized
    return output_acf
end

function anal_acf(output_bit::Array{Int,2}, lags::Int, dt::Float64,
                  p0::Vector{Float64}=[0.7, 3.0, 0.2, 0.1, 0.01])
    # p0 = initial guess
    output_acf = make_acf_dot(output_bit, lags)
    xdata = make_xaxis(lags, dt)
    converged, param, fitting, lifetime = threeExpFit(output_acf, xdata, p0)
    return ACF(converged, param, fitting, lifetime, xdata, [output_acf fitting])
end

"""
Autocorrelation function using dot product
If nr > 1, the binding events for all receptor will be counted for the ACF,
so it will make a single ACF curve. (output: Vector)
"""
function make_acf_dot(output::Array{Int,2}, lags::Int)
    nr, nstep = size(output)
    acf_temp = Array{Float64,2}(undef, (lags + 1, nr))
    temp = Vector{Int}(undef, nstep)
    for i = axes(output, 1)
        temp[:] = view(output, i, :) # copy
        for j = 0:lags
            acf_temp[j + 1, i] = dot(view(temp, 1:(lastindex(temp) - j)),
                                     view(temp, (1 + j):lastindex(temp)))
        end
    end
    acf_sum = sum(acf_temp, dims=2)[:]
    count_num = collect(nstep:-1:(nstep - lags)) * nr
    acf = acf_sum ./ count_num # mean
    acf_norm = acf / acf[1] # normalized
    return acf_norm
end

"""
triple autocorrelation function C3, especially f(x)f(x+τ)f(x+2τ)
If nr > 1, the binding events for all receptor will be counted for the ACF,
so it will make a single ACF curve. (output: Vector)
"""
function make_acf3(output::Array{Int,2}, lags::Int)
    nr, nstep = size(output)
    acf_temp = Array{Float64,2}(undef, (lags + 1, nr))
    temp = Vector{Int}(undef, nstep)
    for i = axes(output, 1)
        temp[:] = view(output, i, :) # copy
        for j = 0:lags
            h1 = view(temp, 1:(lastindex(temp) - 2 * j))
            h2 = view(temp, (1 + j):lastindex(temp) - j)
            h3 = view(temp, (1 + 2 * j):lastindex(temp))
            acf_temp[j + 1, i] = elem_multi_3v(h1, h2, h3)
        end
    end
    acf_sum = sum(acf_temp, dims=2)[:]
    count_num = collect(nstep:-2:(nstep - 2 * lags)) * nr
    acf = acf_sum ./ count_num # mean
    acf_norm = acf / acf[1] # normalized
    return acf_norm
end

""" wating time """
function calc_wait(output::Array{Int,2})
    temp = Vector{Int}(undef, size(output, 2))
    wt_b_u = Vector{Int}() # waiting time from binding to unbinding
    wt_u_b = Vector{Int}() # waiting time from unbinding to binding
    wtime_output = Vector{Float64}(undef, 6)
    for i = axes(output, 1)
        temp[:] = view(output, i, :) # copy
        wt_whole = rle(temp) # run-length encoding of a vector as a tuple
        if wt_whole[1][begin] == 1 # if trajectory starts from binding
            append!(wt_b_u, wt_whole[2][(begin + 2):2:end] .- 1) # remove first
            append!(wt_u_b, wt_whole[2][(begin + 1):2:end] .- 1)
        else # if trajectory starts from unbinding
            append!(wt_b_u, wt_whole[2][(begin + 1):2:end] .- 1)
            append!(wt_u_b, wt_whole[2][(begin + 2):2:end] .- 1) # remove first
        end
        if wt_whole[1][end] == 1 # remove the last one
            pop!(wt_b_u)
        else
            pop!(wt_u_b)
        end
    end
    # wtime_u = fit(Histogram, wt_b_u, minimum(wt_b_u):maximum(wt_b_u) + 1)
    # wtime_b = fit(Histogram, wt_u_b, minimum(wt_u_b):maximum(wt_u_b) + 1)
    m1_b_u = mean(wt_b_u) # first moment (b -> u)
    m2_b_u = mean(wt_b_u.^2) # second moment (b -> u)
    m3_b_u = mean(wt_b_u.^3) # third moment (b -> u)
    m1_u_b = mean(wt_u_b) # first moment (u -> b)
    m2_u_b = mean(wt_u_b.^2) # second moment (u -> b)
    m3_u_b = mean(wt_u_b.^3) # third moment (u -> b)
    wtime_output = [m1_b_u; m2_b_u; m3_b_u; m1_u_b; m2_u_b; m3_u_b]
    return Wtime(wt_b_u, wt_u_b, wtime_output)
end =#


end # module
