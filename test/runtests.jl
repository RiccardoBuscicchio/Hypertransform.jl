using Test
using Hypertransform

@testset "Hypertransform" begin
    x = [0.3, 0.9, 0.1]
    y = hypertriangulate(x)
    x_back = hypercubify(y)
    @test isapprox(x, x_back; atol=1e-10)
    @test size(y) == size(x)

    xmat = [0.1 0.5 0.7; 0.3 0.2 0.6]
    ymat = hypertriangulate(xmat)
    xmat_back = hypercubify(ymat)
    @test isapprox(xmat, xmat_back; atol=1e-10)
    @test size(ymat) == size(xmat)

    # Test sorted output
    @test all(eachrow(ymat)) do row
        issorted(collect(row))
    end

    # Test with a (100, 10) matrix of random numbers
    xrand = rand(100, 10)
    yrand = hypertriangulate(xrand)
    xrand_back = hypercubify(yrand)
    @test isapprox(xrand, xrand_back; atol=1e-10)
    @test size(yrand) == size(xrand)
end