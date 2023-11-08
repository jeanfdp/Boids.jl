using Boids
using Test

using LinearAlgebra: normalize,normalize!,norm

@testset "Boids.jl" begin
    @testset "BoidState function" begin
        n = 10
        box_size = (5, 5)
        boid_state = BoidState(n, box_size)
        @test length(boid_state.positions) == n
        @test length(boid_state.orientations) == n
        @test all([length(pos) == length(box_size) for pos in boid_state.positions])
    end

    @testset "normPeriodic function" begin
        v1 = [1.0, 2.0]
        v2 = [3.0, 4.0]
        v3 = [5.0,4.0]
        box_size = (5, 5)
        @test normPeriodic(v1, v2, box_size) ≈ norm(v1 - v2)
        @test normPeriodic(v1, v3, box_size) ≈ norm(v1 - (v3-[5.0,0.]))
    end

    @testset "normReflect function" begin
        v1 = [1.0, 2.0]
        v2 = [3.0, 4.0]
        v3 = [5.0,4.0]
        @test normReflect(v1, v2) ≈ norm(v1 - v2)
        @test normReflect(v1, v3) ≈ norm(v1 - v3)
    end

    @testset "movePeriodic! function" begin
        box_size = (5, 5)
        boid_state = BoidState([[1.0, 2.0], [3.0, 4.0]], [[1.0, 0.0], [0.0, 1.0]])
        movePeriodic!(boid_state, 1.0, box_size)
        @test boid_state.positions ≈ [[2.0, 2.0], [3.0, 5.0]]
        @test boid_state.orientations ≈ [[1.0, 0.0], [0.0, 1.0]]

        boid_state = BoidState([[1.0, 2.0], [3.0, 4.1]], [[1.0, 0.0], [0.0, 1.0]])
        movePeriodic!(boid_state, 1.0, box_size)
        @test boid_state.positions ≈ [[2.0, 2.0], [3.0, 0.1]]
        @test boid_state.orientations ≈ [[1.0, 0.0], [0.0, 1.0]]
    end

    @testset "moveReflect! function" begin
        box_size = (5, 5)
        boid_state = BoidState([[1.0, 2.0], [3.0, 4.0]], [[1.0, 0.0], [0.0, 1.0]])
        moveReflect!(boid_state, 1.0, box_size)
        @test boid_state.positions ≈ [[2.0, 2.0], [3.0, 5.0]]
        @test boid_state.orientations ≈ [[1.0, 0.0], [0.0, 1.0]]

        boid_state = BoidState([[1.0, 2.0], [3.0, 4.1]], [[1.0, 0.0], [0.0, 1.0]])
        moveReflect!(boid_state, 1.0, box_size)
        @test boid_state.positions ≈ [[2.0, 2.0], [3.0, 4.9]]
        @test boid_state.orientations ≈ [[1.0, 0.0], [0.0, -1.0]]
    end

end