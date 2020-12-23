using Statistics
using LinearAlgebra

function prz(r, z, f, w) # compute poles, residues, zeros
    m = length(w); B = Matrix{Float64}(I, m+1, m+1); B[1, 1] = 0
    E = [0 Transpose(w); ones(m) Diagonal(z)]
    pol = eigvals(E, B); pol = pol[.!isinf.(pol)] # poles
    dz = @. 1e-5*exp(2im*pi*(1:4)/4)
    res = r(pol .+ Transpose(dz))*dz/4 # residues
    E = [0 Transpose(w.*f); ones(m) Diagonal(z)]
    zer = eigvals(E, B); zer = zer[.!isinf.(zer)] # zeros
    return pol, res, zer
end

"""
    r, pol, res, zer, z, f, w, errvec = $(FUNCTIONNAME)(F, Z; tol=1e-13, mmax=100)

AAA rational approximation [^Trefethen2018] of data `F` on set `Z`.

Inputs:
- `F`: vector of data values, or a function handle
- `Z`: vector of sample points
- `tol`: relative tolerance tol
- `mmax`: max type is `(mmax-1,mmax-1)`

Outputs:
- `r`: AAA approximant to F (function)
- `pol,res,zer`: vectors of poles, residues, zeros
- `z,f,w`: vectors of support pts, function values, weights
- `errvec` vector of errors at each step

[^Trefethen2018]: Nakatsukasa, Yuji, Olivier Sète, and Lloyd N. Trefethen. "The AAA algorithm for rational approximation." SIAM Journal on Scientific Computing 40.3 (2018): A1494-A1522.
"""
function aaa(
             F::AbstractVector{T}, Z;
             tol=1e-13,
             mmax=100,
            ) where T
    M = length(Z) # number of sample points
    # `I` is the index set outside the control points
    # `J` is the index set inside the control points
    I = collect(1:M); J = Int[]
    # initializations
    errvec = T[]; w = T[]
    tol′ = norm(F, Inf) * tol

    # select the first support point
    avg = mean(F); err = zero(T); j = 0
    for (i, f) in enumerate(F)
        errᵢ = abs(f - avg)
        if errᵢ > err
            j = i
        end
    end

    # main loop
    for m in 1:mmax
        # update index vector
        push!(J, j)
        deleteat!(I, searchsortedfirst(I, j))

        # solve a linear least-squares problem
        Aᵐ = löwner_matrix(F, Z, I, J)
        _, _, V = svd(Aᵐ, full=true)
        w = V[:, m] # weight vector = min singular value vector

        # max error at sample points
        # i.e. max error at `Z[ii] ∀ ii ∈ I`
        err = zero(T); j = 0
        for ii in I
            Rᵢ = _eval_barycentric(Z[ii], w, F, Z, J)
            errᵢ = abs(F[ii] - Rᵢ)
            if errᵢ > err
                err = errᵢ
                j = ii
            end
        end
        push!(errvec, err)
        err <= tol′ && break # stop if converged
    end
    f, z = F[J], Z[J]
    r = let z=z, f=f, w=w
        zz -> eval_barycentric(zz, w, f, z) # AAA approximant
    end
    pol, res, zer = prz(r, z, f, w) # poles, residues, and zeros
    #r,pol,res,zer,z,f,w = cleanup(r,pol,res,zer,z,f,w,Z,F); # remove Frois. doublets (optional)
    return r, pol, res, zer, z, f, w, errvec
end

function löwner_matrix(F::AbstractVector{T}, Z, I, J) where T
    M = length(Z)
    m̂ = length(I)
    m = length(J)
    @assert m̂ + m == M
    Aᵐ = zeros(T, m̂, m)

    for (j, jj) in enumerate(J), (i, ii) in enumerate(I)
        Aᵐ[i, j] = (F[ii] - F[jj]) / (Z[ii] - Z[jj])
    end
    return Aᵐ
end

# internal use only
# evaluate barycentric form on a subset
function _eval_barycentric(zz, w, f::AbstractVector{T}, z, J) where T
    partial_d = partial_n = zero(T)
    @assert length(w) == length(J)
    for (j, jj) in enumerate(J)
        d = w[j] / (zz - z[jj])
        n = d * f[jj]
        partial_d += d
        partial_n += n
    end
    return partial_n / partial_d
end

eval_barycentric(zz, w, f, z) = map(let w=w, f=f, z=z
                                        zv->eval_barycentric(zv, w, f, z)
                                    end, zz)
function eval_barycentric(zz::Number, w, f, z)
    T = eltype(f)
    partial_d = partial_n = zero(T)
    @assert length(w) == length(f) == length(z)
    @inbounds for j in eachindex(z)
        zz == z[j] && return f[j] # force interpolation
        d = w[j] / (zz - z[j])
        n = d * f[j]
        partial_d += d
        partial_n += n
    end
    return partial_n / partial_d
end
