using NumericalRationalFunctions
using LinearAlgebra
using Test

@testset "AAA" begin
    Z = range(-1, stop=1, length=20)
    F = exp.(Z)
    r, pol, res, zer, z, f, w, errvec = aaa(F, Z)
    @test r(Z) ≈ F atol=1e-14
    @test norm(r(zer), Inf) < 1e-9
    @test errvec[end] < 1e-13
end
