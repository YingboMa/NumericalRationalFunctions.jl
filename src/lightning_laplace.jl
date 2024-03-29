using PolygonOps: inpolygon

function process_domain(
                        P;
                        tol = 1e-6,
                        rel = false,
                        boundary = nothing,
                       )
    randomcirc = false
    if !(P isa AbstractArray && !(eltype(P) <: Number))
        if eltype(P) <: Number
            if length(P) > 1
                w = P # vertices have been specified
            else
                if P < 0
                    randomcirc = true
                    P = -P # random circular arcs
                end
                w = @. exp(2im*pi*(1:P)/P)*(0.1+$rand(P))  # random vertices
            end
        else
            if P == :L
                w = [2.0, 2+1im, 1+1im, 1+2im, 2im, 0.0]
            elseif P == :circleL
                P = complex.([2, [2+1im, -1], 1+2im, 2im, 0.0])
            elseif P == :pent
                w = @. 0.7*exp(pi*2im*(1:5)/5)
            elseif P == :snow
                P = @. exp(2i*pi*(1:12)/12)
                w = @. P*(1+0.2*(-1)^(1:12)) / 1.4
            elseif P == :iso
                w = [1+2im, 1+3im, 2im, 1im+1, 2+1im, 2, 3+1im, 3.0+2im]
                @. w = (w-(1.5+1.5im))/1.8
            end
        end
        !(P isa AbstractArray) && (P = w)
        P = Vector{Any}(P)
        if randomcirc
            for k = 1:length(P)
                r = 0.6/rand()
                P[k] = [P[k]; (rand(Bool) ? r : -r)]
            end
        end
    end

    nw = length(P)
    w = map(first, vec(P))
    ww = similar(w, 0)                                            # bndry pts for plotting
    pt = []
    dw = zeros(nw)
    for k = 1:nw
        kn = mod(k, nw)+1                                # index of next corner
        push!(ww, w[k])
        if eltype(P[k]) <: Number
            if P[k] isa Number                             #     straight arc
                dw[k] = abs(w[kn]-w[k]);                   # distance to next corner
                push!(pt, let w=w, k=k, kn=kn
                          t -> w[k] + t*(w[kn]-w[k])/dw[k]
                      end)  # parametrization of arc
            else                                          #     circular arc
                r = P[k][2];                               # radius of arc
                a = w[k]; b = w[kn]; ab = abs(b-a);        # endpoints of arc
                theta = asin(ab/(2*r));                    # half-angle of arc
                c = a + r*exp(im*(pi/2-theta))*(b-a)/ab;   # center of arc
                dw[k] = 2*theta*r;                         # arc length of arc
                push!(pt, let r=r, theta=theta, a=a, b=b, ab=ab, c=c
                          # parametrization of arc
                          t -> c - r*exp(im*(pi/2+t/r-theta))*(b-a)/ab
                      end)
                append!(ww, pt[k].(range(0,stop=dw[k],length=50)))
            end
        else
            error("general boundary arcs not yet implemented")
        end
    end
    push!(ww, w[1])

    ## Next treat the boundary conditions
    default_bd = z -> real(z).^2
    g = Vector{Any}(undef, nw)
    for k = 1:nw
        g[k] = default_bd       # default
    end
    # This block specifies Dirichlet bndry data g.
    if boundary === nothing
    elseif boundary isa AbstractArray
        if eltype(boundary) <: Function
            g = boundary
        elseif eltype(boundary) <: Number     # if vector, convert to array of fun.
            for k = 1:nw
                g[k] = let boundary = boundary[k]
                    z -> oftype(z, boundary)
                end
            end
        else
            @goto BD_ERROR
        end
    else
        if boundary isa Function
            g = fill(boundary, nw)
        elseif boundary isa Number     # if vector, convert to array of fun.
            for k = 1:nw
                g[k] = let boundary = boundary
                    z -> oftype(z, boundary)
                end
            end
        else
            @goto BD_ERROR
        end
    end

    continuous = true         # check for disc. bndry data if 'rel' not specified
    for k = 1:nw
        j = mod(k-2,nw)+1
        gkk = g[k](w[k])
        gjk = g[j](w[k])
        if abs(gkk-gjk) > tol || isnan(gkk) || isnan(gjk)
            continuous = false   # continuity not enforced at Neumann corners
        end
    end
    if !continuous
        rel = true
    end

    return g, w, ww, pt, dw, tol, rel

    @label BD_ERROR
    error("Invalid boundary. Only use functions or numbers.")
end

function laplace(
                 P;
                 tol = 1e-6,
                 steps = false,
                 plots = true,
                 slow = false,
                 rel = false,
                 aaa = true,
                 boundary = nothing,
                 arnoldi = true,
                )
    steps && (plots = true)
    g, w, ww, pt, dw, tol, rel = process_domain(P; tol, rel, boundary)
    Zplot = ww

    # number of corners
    nw = length(w)
    wr = extrema(real(ww))
    wi = extrema(imag(ww))
    # for scale- and transl-invariance
    wc = mean(@. wr + wi * im)
    scl = max(wr[2] - wr[1], wi[2] - wi[1])
    # sets which corners get more poles
    q = 0.5
    slow && (q = 0.0)

    inpolygonc = let polygon=map(w->(real(w), imag(w)), ww)
        zs -> map(z -> inpolygon((real(z), imag(z)), polygon), zs)            # complex variant of "inpolygon"
    end

    TT = complex(float(eltype(ww)))
    outward = similar(ww, TT)
    for k = 1:nw
        forward = pt[k](0.01*dw[k]) - w[k]            # small step toward next corner
        j = mod(k-2,nw)+1
        backward = pt[j](0.99*dw[j]) - w[k]           # small step toward last corner
        tmp = im*backward*sqrt(-forward/backward);
        outward[k] = tmp/abs(tmp);                    # outward direction from corner
    end

    Nvec = Int[]
    errvec = real(TT)[]
    errk = ones(nw)
    nkv = zeros(Int, nw)
    maxstepno = 30
    err0 = Inf
    err = real(TT)

    Z = TT[]            # vector of sample points on boundary
    G = TT[]            # vector of boundary values at these points
    T = Vector{TT}[]    # vector of unit tangent vectors at these points
    pol = TT[]          # vector of poles of the rational approximation
    J = Int[]           # vector of indices of which corner each pole belongs to
    d = real(TT)[]      # vector of distances from poles to their corners
    tt = map(_->real(TT)[], 1:nw)  # array of distances of sample points along each side

    # initializations
    cc = pol
    A = Matrix{real(TT)}(undef, 0, 0)
    H = Matrix{TT}(undef, 0, 0)
    n = 0
    for stepno = 1:maxstepno
        empty!(Z); empty!(G); empty!(T); empty!(pol); empty!(J); empty!(d)
        for t in tt
            empty!(t)
        end
        # Fix poles and sample pts on bndry.  Side k means side from corner k to k+1.
        for k = 1:nw
            nk = nkv[k]                                  # no. of poles at this corner
            sk = @. sqrt(1:nk) - $sqrt(nk)
            dk = @. exp(4*sk); dk = scl*dk;                  # stronger clustering near corner
            dk = dk[@. dk > 1e-15*scl]                        # remove poles too close to corner
            polk = @. w[k] + outward[k]*dk;                  # poles near this corner
            ii = findfirst(isequal(1), inpolygonc(polk[@. dk>1e-12*scl])) # work around inaccuracy
            if ii !== nothing                               # don't allow poles in Omega
                dk = dk[1:ii-2]
                polk = polk[1:ii-2]
            end
            append!(pol, polk); append!(d, dk)
            for _ in eachindex(dk)
                push!(J, k)
            end
            dvec = [(1/3)*dk; (2/3)*dk; dk]                # finer pts for bndry sampling
            clustered = dvec[dvec.<dw[k]] # clustered pts near corner
            additional = range(0,stop=dw[k],length=max(30,nk))  # additional pts along side
            append!(tt[k], clustered)
            append!(tt[k], additional)
            j = mod(k-2,nw)+1                            # index of last corner
            append!(tt[j], dw[j].-dvec[dvec.<dw[j]])       # likewise in other direction
        end
        for k = 1:nw
            sort!(tt[k])
            tk = tt[k]; pk = pt[k]                       # abbrevations
            pktk = pk.(tk)
            append!(Z, pktk)                             # sample pts on side k
            append!(G, g[k].(pktk))                        # boundary data at these pts
            h = 1e-4                                     # 4-pt trapezoidal rule
            push!(T, @. (pk.(tk.+h)-im*pk.(tk.+im*h) - pk(tk-h)+im*pk(tk-im*h))/(4*h))    # unnormalized tangent vectors
        end
        for t in T
            @. t = t/abs(t)                                  # normalize tangent vectors
        end
        II = isnan.(G);                                   # Neumann indices
        if any(II)
            arnoldi = false
        end

        # Solve the Laplace problem
        n = 4stepno                                    # degree of polynomial term
        Np = length(pol)

        M = size(Z, 1)
        H = zeros(TT, n+1, n);                                # Arnoldi Hessenberg matrix
        if arnoldi
            Q = ones(TT, M, n+1) # cols of Q have norm sqrt(M)
            v = zeros(TT, M)
            @views for k = 1:n                                   # Arnoldi process
                @. v = (Z - wc) * Q[:, k]                        # (not yet implemented if
                for j = 1:k                                #  are Neumann BCs)
                    H[j, k] = Q[:, j]'v/M
                    @. v = v - H[j, k]*Q[:, j]
                end
                H[k+1, k] = norm(v)/sqrt(M)
                @. Q[:, k+1] = v/H[k+1, k]
            end
        else                                             # no-Arnoldi option
            Q = @. ((Z-wc)/scl)^$Transpose(0:n)                      # (for educational purposes)
        end
        A = [real(Q) imag(Q[:,2:n+1])]                  # matrix for least-sq
        if any(II)                                       # Neumann BCs, if any
            tmp = @. (0:n)*((Z[II]-wc)/scl).^[0; 0:n-1]*T[II]/scl
            A[II,:] = [imag(tmp) -real(tmp[:,2:n+1])]
        end
        # Dirichlet row entry pairs have the form Re (u+iv)*(a-ib) = [u v][a b]' = g
        # Neumann row entry pairs have the form Im (U+iV)*(a-ib) = [V -U][a b]' = 0
        # where U+iV = (u+iv)/(unit vector T in tangent direction)
        # and now u and v correspond not to d/(Z-pol) but to -d/(Z-pol)^2
        if Np > 0                                         # g linear => no poles
            complexA = @. $Transpose(d)/(Z-$Transpose(pol))
            updateA = [real(complexA) imag(complexA)]
            if isempty(A)
                A = updateA
            else
                A = [A updateA]     # columns with poles
            end
            if any(II)                                     # Neumann BCs, if any
                JJ = 2*n+1 .+ (1:2*Np);                      # column indices for poles
                A[II,1] = 0;
                tmp = @. -(d/(Z[II]-pol)^2)*T[II]
                A[II,JJ] = [imag(tmp) -real(tmp)]
            end
        end
        # corner for each col
        for _ in 1:2n+1
            push!(J, 0)
        end
        append!(J, J); append!(J, J)
        N = size(A,2)                                    # no. of cols = 2n+1+2Np
        Kj = zeros(Int, M)
        for j = 1:M
            dd = @. abs(Z[j]-w)
            Kj[j] = findfirst(isequal(minimum(dd)), dd)                   # nearest corner to Zj
        end
        if rel                                            # weights to measure error
            wt = @. abs(Z-w[Kj])/scl
        else
            wt = ones(M)
        end
        W = Diagonal(sqrt.(wt))                      # weighting for case 'rel'
        Gn = G; Gn[II] .= 0                               # set Neumann vals to 0
        c = (W*A)\(W*Gn)                                 # least-squares solution
        cc = [c[1]; c[2:n+1]-im*c[n+2:2*n+1]              # complex coeffs for f
              c[2*n+2:2*n+Np+1]-im*c[2*n+Np+2:end]]
        for k = 1:nw
            Kk = findall(isequal(k), Kj)
            errk[k] = norm(wt[Kk] .* (A[Kk,:]*c-Gn[Kk]), Inf) # error near corner k
        end
        err = norm(wt.*(A*c-Gn), Inf)                     # global error
        polmax = 100
        for k = 1:nw
            if (errk[k] > q*err) & (nkv[k] < polmax)
                nkv[k] = nkv[k]+ceil(1+sqrt(nkv[k]))      # increase no. poles
            else
                nkv[k] = max(nkv[k],ceil(stepno/2))
            end
        end
        push!(errvec, err); push!(Nvec, N)
        err < 0.5*tol && break                             # convergence success
        if err < err0                                      # save the best so far
            #TODO
        end
        if (N > 1200) || (stepno == maxstepno) || (Np >= polmax*nw)  # failure
            #TODO
            #u = u0; f = f0; Z = Z0; G = G0; A = A0; M = M0;
            #N = N0; err = err0; pol = pol0; wt = wt0;
            @warn "Loosen tolerance or add corners?"
            break
        end
    end
    f = let wc=wc, cc=cc, H=H, pol=pol, d=d, arnoldi=arnoldi, scl=scl, n=n
        z -> fzeval(z,wc,               # vector and matrix inputs
                    cc,H,pol,d,arnoldi,scl,n)    # to u and f both allowed
    end
    u = let f=f
        z -> real(f(z))
    end
    return u, err, f, Z, Zplot, A, inpolygonc
end

function fzeval(Z,wc,cc,H,pol,d,arnoldi,scl,n)
    !(Z isa Number) && (Z = vec(Z))
    ZZ = [wc; Z]
    if arnoldi
        M = length(ZZ)
        TT = eltype(H)
        Q = ones(TT, M, n+1)
        v = zeros(TT, M)
        @views for k = 1:size(H,2)
            @. v = (ZZ-wc)*Q[:,k]
            for j = 1:k
                @. v = v - H[j,k]*Q[:,j]
            end
            @. Q[:, k+1] = v/H[k+1, k]
        end
    else
        Q = @. ((ZZ-wc)/scl).^$Transpose(0:n);
    end
    if length(pol) > 0
        fZZ = [Q (@. $Transpose(d)/(ZZ-$Transpose(pol)))]*cc
    else
        fZZ = Q*cc
    end
    fZ = @. fZZ[2:end] - im*imag(fZZ[1])
    return Z isa Number ? fZ[1] : reshape(fZ, size(Z))
end
