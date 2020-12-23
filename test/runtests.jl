using NumericalRationalFunctions
using LinearAlgebra
using Test

@testset "AAA" begin
    Z = range(-1, stop=1, length=20)
    F = exp.(Z)
    r, pol, res, zer, z, f, w, errvec = aaa(F, Z)
    @test length(z) == 7
    @test length(pol) == length(res) == length(zer) == 6
    @test r(Z) ≈ F atol=1e-14
    x = 0.12345
    @test r(x) ≈ exp(x) atol=1e-14
    @test norm(r(zer), Inf) < 1e-9
    @test errvec[end] < 1e-13
end

@testset "Laplace" begin
    u, err, f, Z, Zplot, A, inpoly = laplace(:iso, tol=1e-10, boundary=z->log(abs(z)))
    @test err < 1e-10
    u, err, f, Z, Zplot, A, inpoly = laplace(:circleL)
    @test err < 1e-6
    u, err, f, Z, Zplot, A, inpoly = laplace(complex.([[1, 0.5], 1+1im, -1+1im, -1]), boundary=[1, 0, 0, 0])
    @test err < 1e-6

    #=
    sx = sy = range(-1.5, stop=1.5, length=200)
    p = contour(sx, sy, (x, y)->begin
                    z = x + im*y
                    inpoly(z) == 1 ? u(z) : NaN
                end, levels=30, c=:viridis, aspect_ratio=1)
    nw = length(Zplot)
    for k in 1:nw-1
        kn = k+1
        plot!(p, [Zplot[k], Zplot[kn]], color=:black, lab=false, linewidth=2)
    end
    p
    =#
end
