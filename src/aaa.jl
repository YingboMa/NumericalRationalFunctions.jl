using Statistics
using LinearAlgebra

function prz(r, z, f, w) # compute poles, residues, zeros
    m = length(w); B = Matrix{Float64}(I, m+1, m+1); B[1, 1] = 0
    E = [0 Transpose(w); ones(m) Diagonal(z)];
    pol = eigvals(E, B); pol = pol[.!isinf.(pol)] # poles
    dz = @. 1e-5*exp(2im*pi*(1:4)/4);
    res = r(pol .+ Transpose(dz))*dz/4; # residues
    E = [0 Transpose(w.*f); ones(m) Diagonal(z)];
    zer = eigvals(E, B); zer = zer[.!isinf.(zer)]; # zeros
    return pol, res, zer
end

function rhandle(zz, z, f, w) # evaluate r at zz
    zv = zz isa Number ? zz : vec(zz)
    CC = 1 ./ (zv .- Transpose(z)); # Cauchy matrix
    r = (CC*(w.*f))./(CC*w) # AAA approx as vector
    ii = findall(isnan.(r)) # find values NaN = Inf/Inf if any
    for j = 1:length(ii)
        r[ii[j]] = f[findfirst(isequal(zv[ii[j]]), z)] # force interpolation there
    end
    # AAA approx.
    if zz isa Number
        return r[1]
    else
        return reshape(r, size(zz))
    end
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

!!! note
    The code is a direct translation of the code presented in the paper with
    minor optimization.
"""
function aaa(
             F::AbstractVector{T}, Z;
             tol=1e-13,
             mmax=100,
            ) where T
    M = length(Z) # number of sample points
    SF = Diagonal(F) # left scaling matrix
    J = collect(1:M); z = T[]; f = T[]; C = T[]; # initializations
    errvec = T[]; R = mean(F);
    w = T[]
    for m = 1:mmax # main loop
        j = argmax(@. abs(F-R)) # select next support point
        push!(z, Z[j])
        f = [f; F[j]] # update support points, data values
        deleteat!(J, J .== j) # update index vector
        if isempty(C)
            C = reshape((@. 1/(Z-Z[j])), M, 1) # next column of Cauchy matrix
        else
            C = [C (@. 1/(Z-Z[j]))] # next column of Cauchy matrix
        end
        Sf = Diagonal(f) # right scaling matrix
        A = SF*C - C*Sf # Loewner matrix
        _, _, V = svd(A[J,:], full=true) # SVD
        w = V[:, m] # weight vector = min sing vector
        N = C*(w.*f); D = C*w; # numerator and denominator
        R = copy(F); @. R[J] = N[J]/D[J] # rational approximation
        err = norm(F-R, Inf)
        push!(errvec, err) # max error at sample points
        err <= tol*norm(F, Inf) && break # stop if converged
    end
    r = let z=z, f=f, w=w
        zz -> rhandle(zz, z, f, w) # AAA approximant as function handle
    end
    pol, res, zer = prz(r, z, f, w) # poles, residues, and zeros
    #r,pol,res,zer,z,f,w = cleanup(r,pol,res,zer,z,f,w,Z,F); # remove Frois. doublets (optional)
    return r, pol, res, zer, z, f, w, errvec
end
