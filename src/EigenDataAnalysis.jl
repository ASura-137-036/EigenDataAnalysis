module EigenDataAnalysis

using LinearAlgebra
using Statistics

export SystemIDResult, SSIResult, StabilizationPoint
export SSICOV, SSIDATA, PCSSI, SystemIDMethod
export fit, mac, filter_modes, reconstruct, damping_ratios, compute_stabilization_diagram, reconstruction_error
export singular_values, system_matrices, modes, eigenvalues

# --- Abstract Types & Interfaces ---
abstract type SystemIDMethod end
struct SSICOV <: SystemIDMethod end
struct SSIDATA <: SystemIDMethod end
struct PCSSI <: SystemIDMethod end

abstract type SystemIDResult end

# --- Unified Result Container ---
struct SSIResult{T<:SystemIDMethod} <: SystemIDResult
    method_type::Type{T}
    dt::Float64
    continuous_eigenvalues::Vector{ComplexF64}
    modes::Matrix{ComplexF64}
    amplitudes::Vector{ComplexF64}
    A::Matrix{Float64}
    C::Matrix{Float64}
    S::Vector{Float64}
    max_lag::Int
    state_sequence::Matrix{Float64}   # Size: (r_stable, N_cols) -> Empty for SSICOV
    output_sequence::Matrix{Float64}  # Size: (ny, N_cols) or (ny, N)
    time_indices::UnitRange{Int}      # Exact time indices for proper alignment
end

# --- Accessor API ---
singular_values(res::SSIResult) = res.S
system_matrices(res::SSIResult) = (A = res.A, C = res.C)
modes(res::SSIResult) = res.modes
eigenvalues(res::SSIResult) = res.continuous_eigenvalues

# --- Analytical Amplitude Solver (SSICOV) ---
function compute_global_amplitudes(Φ::AbstractMatrix, ω::Vector{ComplexF64}, X::AbstractMatrix, dt::Float64)
    ny, N = size(X)
    r = length(ω)
    if r == 0
        return ComplexF64[]
    end
    t = (0:N-1) .* dt
    
    M = zeros(ComplexF64, ny * N, r)
    for k in 1:r
        # Uses transpose() to avoid complex conjugation of frequencies
        spatial_temporal = Φ[:, k] * transpose(exp.(ω[k] .* t))
        M[:, k] = vec(spatial_temporal)
    end
    return M \ vec(ComplexF64.(X))
end

# ==============================================================================
# METHOD 1: SSI-COV (Block-Hankel of Covariances)
# ==============================================================================
function fit(::Type{SSICOV}, X::AbstractMatrix, dt::Float64; max_lag::Int=20, r::Int=10 , unstable_filt::Bool=false)
    ny, N = size(X)
    
    R = [zeros(ny, ny) for _ in 1:(2*max_lag-1)]
    for k in 1:(2*max_lag-1)
        R[k] = (X[:, 1+k:N] * X[:, 1:N-k]') ./ (N - k)
    end
    
    H = zeros(ny * max_lag, ny * max_lag)
    for i in 1:max_lag, j in 1:max_lag
        H[(i-1)*ny+1 : i*ny, (j-1)*ny+1 : j*ny] = R[i + j - 1]
    end
    
    U, S, V = svd(H)
    r_actual = min(r, length(S))
    O = U[:, 1:r_actual] * Diagonal(sqrt.(S[1:r_actual]))
    
    C_mat = O[1:ny, :]
    A_mat = O[1:end-ny, :] \ O[ny+1:end, :]
    
    Λ, Ψ = eigen(A_mat)
    Φ_full = C_mat * Ψ
    ω_full = log.(Complex.(Λ)) ./ dt
    
    if unstable_filt == true
        # Filter out physically unstable poles (Re(ω) > 0)
        stable_idx = findall(x -> real(x) <= 0.0, ω_full)
    else
        stable_idx = collect(1:length(ω_full))
        #warn about unstable poles if any exist , with all info
        if !isempty(findall(x -> real(x) > 0.0, ω_full))
            @warn "Unstable poles detected in SSICOV fit. Consider using unstable_filt=true to filter them out."
            for idx in findall(x -> real(x) > 0.0, ω_full)
                @info "Unstable pole: ω=$(ω_full[idx]), |Λ|=$(abs(Λ[idx])), mode=$(Φ_full[:, idx])"
            end
        end
    end
    ω = ω_full[stable_idx]
    Φ = Φ_full[:, stable_idx]
    
    b = compute_global_amplitudes(Φ, ω, X, dt)
    
    t = (0:N-1) .* dt
    output_seq = isempty(ω) ? zeros(ny, N) : real.(Φ * Diagonal(b) * exp.(ω * transpose(t)))
    
    return SSIResult(SSICOV, dt, ω, Φ, b, A_mat, C_mat, S, max_lag, 
                     zeros(0,0), output_seq, 1:N)
end

# ==============================================================================
# METHOD 2: PC-SSI (Block-Toeplitz of Past/Future Data)
# ==============================================================================
function fit(::Type{PCSSI}, X::AbstractMatrix, dt::Float64; max_lag::Int=20, r::Int=10 , unstable_filt::Bool=false)
    ny, N = size(X)
    N_cols = N - 2*max_lag + 1
    
    Yp = zeros(ny * max_lag, N_cols)
    Yf = zeros(ny * max_lag, N_cols)
    
    for i in 1:max_lag
        Yp[(i-1)*ny+1 : i*ny, :] = X[:, i : N_cols + i - 1]
        Yf[(i-1)*ny+1 : i*ny, :] = X[:, i + max_lag : N_cols + i + max_lag - 1]
    end
    
    T_mat = (Yf * Yp') ./ N_cols
    U, S, _ = svd(T_mat)
    r_actual = min(r, length(S))
    S_r = S[1:r_actual]
    
    O = U[:, 1:r_actual] * Diagonal(sqrt.(S_r))
    C_mat = O[1:ny, :]
    
    X_seq_actual = O \ Yf
    
    X_k   = X_seq_actual[:, 1:end-1]
    X_kp1 = X_seq_actual[:, 2:end]
    A_mat = X_kp1 / X_k
    
    Λ, Ψ = eigen(A_mat)
    Φ_full = C_mat * Ψ
    ω_full = log.(Complex.(Λ)) ./ dt
    
    if unstable_filt == true
        # Filter out physically unstable poles (Re(ω) > 0)
        stable_idx = findall(x -> real(x) <= 0.0, ω_full)
    else
        stable_idx = collect(1:length(ω_full))
        #warn about unstable poles if any exist , with all info
        if !isempty(findall(x -> real(x) > 0.0, ω_full))
            @warn "Unstable poles detected in PC-SSI fit. Consider using unstable_filt=true to filter them out."
            for idx in findall(x -> real(x) > 0.0, ω_full)
                @info "Unstable pole: ω=$(ω_full[idx]), |Λ|=$(abs(Λ[idx])), mode=$(Φ_full[:, idx])"
            end
        end
    end
    ω = ω_full[stable_idx]
    Φ = Φ_full[:, stable_idx]
    
    # Project state sequence onto stable modal coordinates
    Z_seq_full = Ψ \ X_seq_actual
    Z_seq_stable = Z_seq_full[stable_idx, :]
    
    b = isempty(stable_idx) ? ComplexF64[] : Z_seq_stable[:, 1]
    output_seq = isempty(stable_idx) ? zeros(ny, N_cols) : real.(Φ * Z_seq_stable)
    X_seq_stable = isempty(stable_idx) ? zeros(0, N_cols) : real.(Ψ[:, stable_idx] * Z_seq_stable)
    
    t_indices = (max_lag + 1):(max_lag + N_cols)
    
    return SSIResult(PCSSI, dt, ω, Φ, b, A_mat, C_mat, S, max_lag, 
                     X_seq_stable, output_seq, t_indices)
end

# ==============================================================================
# METHOD 3: SSI-DATA (Orthogonal Projection via LQ Decomposition)
# ==============================================================================
function fit(::Type{SSIDATA}, X::AbstractMatrix, dt::Float64; max_lag::Int=20, r::Int=10 , unstable_filt::Bool=false)
    ny, N = size(X)
    N_cols = N - 2*max_lag + 1
    
    Yp = zeros(ny * max_lag, N_cols)
    Yf = zeros(ny * max_lag, N_cols)
    
    for i in 1:max_lag
        Yp[(i-1)*ny+1 : i*ny, :] = X[:, i : N_cols + i - 1]
        Yf[(i-1)*ny+1 : i*ny, :] = X[:, i + max_lag : N_cols + i + max_lag - 1]
    end
    
    F = lq([Yp; Yf])
    L = F.L
    Q = F.Q
    
    ny_block = ny * max_lag
    L21 = L[ny_block+1 : end, 1:ny_block]
    Q1 = Q[1:ny_block, :]
    
    U, S, V = svd(L21)
    r_actual = min(r, length(S))
    S_r = S[1:r_actual]
    
    O = U[:, 1:r_actual] * Diagonal(sqrt.(S_r))
    C_mat = O[1:ny, :]
    A_mat = O[1:end-ny, :] \ O[ny+1:end, :]
    
    X_seq_actual = Diagonal(sqrt.(S_r)) * V[:, 1:r_actual]' * Q1
    
    Λ, Ψ = eigen(A_mat)
    Φ_full = C_mat * Ψ
    ω_full = log.(Complex.(Λ)) ./ dt
    
    if unstable_filt == true
        # Filter out physically unstable poles (Re(ω) > 0)
        stable_idx = findall(x -> real(x) <= 0.0, ω_full)
    else
        stable_idx = collect(1:length(ω_full))
        #warn about unstable poles if any exist , with all info
        if !isempty(findall(x -> real(x) > 0.0, ω_full
))
            @warn "Unstable poles detected in SSIDATA fit. Consider using unstable_filt=true to filter them out."
            for idx in findall(x -> real(x) > 0.0, ω_full)
                @info "Unstable pole: ω=$(ω_full[idx]), |Λ|=$(abs(Λ[idx])), mode=$(Φ_full[:, idx])"
            end
        end
    end
    ω = ω_full[stable_idx]
    Φ = Φ_full[:, stable_idx]
    
    # Project state sequence onto stable modal coordinates
    Z_seq_full = Ψ \ X_seq_actual
    Z_seq_stable = Z_seq_full[stable_idx, :]
    
    b = isempty(stable_idx) ? ComplexF64[] : Z_seq_stable[:, 1]
    output_seq = isempty(stable_idx) ? zeros(ny, N_cols) : real.(Φ * Z_seq_stable)
    X_seq_stable = isempty(stable_idx) ? zeros(0, N_cols) : real.(Ψ[:, stable_idx] * Z_seq_stable)
    
    t_indices = (max_lag + 1):(max_lag + N_cols)
    
    return SSIResult(SSIDATA, dt, ω, Φ, b, A_mat, C_mat, S, max_lag, 
                     X_seq_stable, output_seq, t_indices)
end

# --- Output Reconstruction & Error ---
function reconstruct(model::SSIResult)
    return model.output_sequence, model.time_indices
end

function reconstruction_error(model::SSIResult, X_original::AbstractMatrix)
    X_align = X_original[:, model.time_indices]
    if size(X_align) != size(model.output_sequence)
        error("Dimension mismatch: Data is $(size(X_align)), model output is $(size(model.output_sequence)).")
    end
    return norm(X_align - model.output_sequence) / norm(X_align)
end

# --- Validation and Diagram Tracking ---
function mac(ϕ1::AbstractVector, ϕ2::AbstractVector)
    num = abs(ϕ1' * ϕ2)^2
    den = real((ϕ1' * ϕ1) * (ϕ2' * ϕ2))
    return real(num / den)
end

function damping_ratios(model::SSIResult)
    ω = model.continuous_eigenvalues
    return -real.(ω) ./ abs.(ω)
end

struct StabilizationPoint
    order::Int
    freq_hz::Float64
    zeta::Float64
    energy::Float64
    status::Symbol
end

function compute_stabilization_diagram(
    Method::Type{<:SystemIDMethod}, X::AbstractMatrix, dt::Real;
    orders::AbstractVector{Int}, tol_f::Float64=0.01, tol_zeta::Float64=0.05, tol_mac::Float64=0.95, kwargs...
)
    points = StabilizationPoint[]
    recon_errors = Float64[]
    prev_freqs, prev_zetas, prev_modes = Float64[], Float64[], Matrix{ComplexF64}(undef, 0, 0)
    
    for r in orders
        model = fit(Method, X, Float64(dt); r=r, kwargs...)
        push!(recon_errors, reconstruction_error(model, X))
        
        ω = eigenvalues(model)
        Φ = modes(model)
        amps = model.amplitudes
        
        valid_idx = findall(x -> imag(x) > 0, ω) # Only positive imaginary frequencies
        curr_freqs = imag.(ω[valid_idx]) ./ (2π)
        curr_zetas = -real.(ω[valid_idx]) ./ abs.(ω[valid_idx])
        curr_modes = Φ[:, valid_idx]
        curr_energies = abs.(amps[valid_idx])
        
        for i in 1:length(valid_idx)
            f_i, z_i, e_i, phi_i = curr_freqs[i], curr_zetas[i], curr_energies[i], curr_modes[:, i]
            status = :new
            
            if !isempty(prev_freqs)
                f_diffs = abs.(prev_freqs .- f_i)
                best = argmin(f_diffs)
                
                is_f_stable = (f_diffs[best] / max(f_i, 1e-5)) < tol_f
                is_z_stable = abs(prev_zetas[best] - z_i) / max(z_i, 1e-5) < tol_zeta
                is_mac_stable = mac(phi_i, prev_modes[:, best]) > tol_mac
                
                if is_f_stable && is_z_stable && is_mac_stable
                    status = :fully_stable
                elseif is_f_stable && is_z_stable
                    status = :freq_damp
                elseif is_f_stable && is_mac_stable
                    status = :freq_mac
                elseif is_f_stable
                    status = :freq
                else
                    status = :unstable
                end
            end
            push!(points, StabilizationPoint(r, f_i, z_i, e_i, status))
        end
        prev_freqs, prev_zetas, prev_modes = curr_freqs, curr_zetas, curr_modes
    end
    return points, recon_errors, orders
end

end # module