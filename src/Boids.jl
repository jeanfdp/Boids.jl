module Boids

using LinearAlgebra: normalize,normalize!,norm
using Distributions: MvNormal

export BoidState, BoidState, update!, ruleOrientations!, rulePosition!, moveReflect!, movePeriodic!, normReflect, normPeriodic, runSimple2DBoids, runSimple3DBoids

"""
    BoidState(positions::Vector{Vector{Float64}}, orientations::Vector{Vector{Float64}})

A structure representing the state of a flock of boids.

# Fields
- `positions`: A vector of vectors representing the positions of each boid in the flock.
- `orientations`: A vector of vectors representing the orientations of each boid in the flock.

# Example
```julia
BoidState([[1.0, 2.0], [3.0, 4.0]], [[0.5, 0.5], [0.7, 0.3]])
```
"""
struct BoidState
    positions::Vector{Vector{Float64}}
    orientations::Vector{Vector{Float64}}
end

"""
    BoidState(n::Int, box_size::Tuple)

Create a new `BoidState` with `n` boids in a box of size `box_size`.

The positions of the boids are randomly initialized within the box, and the orientations are randomly initialized to point in a uniformly random direction.

# Arguments
- `n`: The number of boids.
- `box_size`: A tuple representing the size of the box in each dimension.

# Returns
- A new `BoidState` object.

# Example
In 2D
```julia
BoidState(10, (5, 5))
```

In 3D
```julia
BoidState(10, (5, 5, 5))
```
"""
function BoidState(n::Int, box_size::Tuple)
    dim = length(box_size)
    positions = [rand(dim) .* box_size for _ in 1:n]
    orientations = [normalize(rand(MvNormal(dim, 1))) for _ in 1:n]
    BoidState(positions, orientations)
end

"""
    update!(state::BoidState, ruleOrientations!::Function, rulePositions!::Function, move!::Function)

Update the state of the boids by applying the given rules and movement function.

# Arguments
- `state`: The current state of the boids.
- `ruleOrientations!`: A function that takes only a BoidState and updates the orientations of the boids.
- `rulePositions!`: A function that takes only a BoidState and updates the positions of the boids.
- `move!`: A function that takes only a BoidState and moves the boids according to their updated positions and orientations.

# Example
```julia
update!(newstate,state->ruleOrientations!(state,paramMaxRad,paramVelWeight,(v1,v2)->normPeriodic(v1,v2,(size,size,size))),state->rulePosition!(state,paramMaxRad,paramMinRad,paramPosWeight,paramRepWeight,(v1,v2)->normPeriodic(v1,v2,(size,size,size))),state->movePeriodic!(state,vel,(size,size,size)))
```
"""
function update!(state::BoidState,ruleOrientations!::Function,rulePositions!::Function,move!::Function)
    ruleOrientations!(state)
    rulePositions!(state)
    move!(state)
end


"""
    ruleOrientations!(boid_state::BoidState, paramMaxRad::Float64, paramVelWeight::Float64, norm::Function)

A specific implementation of an alignment rule. Updates the orientations of the boids in `boid_state` based on their neighbors within a radius of `paramMaxRad`.

Each boid's new orientation is shifter by the average of the orientations of its neighbors, where the weight of the shift is given by `paramVelWeight`.

# Arguments
- `boid_state`: The current state of the boids.
- `paramMaxRad`: The maximum radius within which to consider other boids as neighbors.
- `paramVelWeight`: The weight dictating how much a boid's orientation is affected by its neighbors' orientations.
- `norm`: A function to compute the distance between two boids. This function should take two vectors representing the positions of the boids and return a scalar. This is intended to allow both periodic and reflecting norms.

# Example
```julia
ruleOrientations!(boid_state, 1.0, 0.5, normReflect)
```
"""
function ruleOrientations!(boid_state::BoidState, paramMaxRad::Float64, paramVelWeight::Float64,norm::Function)
    for i in eachindex(boid_state.orientations)
        temp = zeros(length(boid_state.orientations[i]))
        flag=false
        for j in eachindex(boid_state.orientations)
            if i == j
                continue
            end
            if norm(boid_state.positions[i], boid_state.positions[j]) > paramMaxRad
                continue
            end
            temp .+= boid_state.orientations[j]
            flag=true
        end
        if !flag
            continue
        end
        normalize!(temp)
        temp .-= boid_state.orientations[i]
        boid_state.orientations[i] .+= temp * paramVelWeight
        normalize!(boid_state.orientations[i])
    end
end

"""
    rulePosition!(boid_state::BoidState, paramMaxRad::Float64, paramMinRad::Float64, paramPosWeight::Float64, paramRepWeight::Float64, norm::Function)

A specific implementation of cohesion and repulsion rules. Updates the positions of the boids in `boid_state` based on their neighbors within a radius of `paramMaxRad`.

Each boid's new position is a weighted average of the positions of its neighbors, where the weight is given by `paramPosWeight`. Boids within a radius of `paramMinRad` are considered too close and are repelled with weight `paramRepWeight`.

# Arguments
- `boid_state`: The current state of the boids.
- `paramMaxRad`: The maximum radius within which to consider other boids as neighbors.
- `paramMinRad`: The minimum radius within which other boids are considered too close and are repelled.
- `paramPosWeight`: The weight to give the effect of the average position of the neighbors when changing direction towards it.
- `paramRepWeight`: The weight to give to the repulsion from boids that are too close.
- `norm`: A function to compute the distance between two boids.

# Example
```julia
rulePosition!(boid_state, 1.0, 0.5, 0.5, 0.5, normPeriodic)
```
"""
function rulePosition!(boid_state::BoidState, paramMaxRad::Float64, paramMinRad::Float64, paramPosWeight::Float64, paramRepWeight::Float64, norm::Function)
    for i in eachindex(boid_state.orientations)
        temp=zeros(length(boid_state.orientations[i]))
        N = 0
        temprepulse = zeros(length(boid_state.orientations[i]))
        Nr = 0
        for j in eachindex(boid_state.positions)
            if i == j
                continue
            end
            if norm(boid_state.positions[i], boid_state.positions[j]) > paramMaxRad
                continue
            end
            if norm(boid_state.positions[i], boid_state.positions[j]) < paramMinRad
                temprepulse .+= boid_state.positions[j]#TODO: This is wrong for periodic boundary conditions

                Nr += 1
                continue
            end
            temp .+= boid_state.positions[j]
            N += 1
        end
        if N != 0
            temp /= N
            temp .-= boid_state.positions[i]
            normalize!(temp)
            boid_state.orientations[i] .+= temp * paramPosWeight
            normalize!(boid_state.orientations[i])
        end
        if Nr != 0
            temprepulse /= Nr
            temprepulse .-= boid_state.positions[i]
            normalize!(temprepulse)
            boid_state.orientations[i] .-= temprepulse * paramRepWeight
            normalize!(boid_state.orientations[i])
        end
    end
end

"""
    moveReflect!(boid_state::BoidState, paramSpeed::Float64, box_size::Tuple)

Move the boids in `boid_state` according to their orientations and `paramSpeed`. If a boid hits the edge of the box, it is reflected back into the box and its orientation is reversed.

# Arguments
- `boid_state`: The current state of the boids.
- `paramSpeed`: The speed at which the boids move.
- `box_size`: A tuple representing the size of the box in each dimension.

# Example
```julia
moveReflect!(boid_state, 1.0, (5, 5))
```
"""
function moveReflect!(boid_state::BoidState, paramSpeed::Float64, box_size::Tuple)
    for i in eachindex(boid_state.positions)
        boid_state.positions[i] .+= boid_state.orientations[i] * paramSpeed
        for j in eachindex(boid_state.positions[i])
            if boid_state.positions[i][j] > box_size[j]
                boid_state.positions[i][j] = 2*box_size[j] - boid_state.positions[i][j]
                boid_state.orientations[i][j] *= -1
            end
            if boid_state.positions[i][j] < 0
                boid_state.positions[i][j] *= -1
                boid_state.orientations[i][j] *= -1
            end
        end
    end
end

"""
    movePeriodic!(boid_state::BoidState, paramSpeed::Float64, box_size::Tuple)

Move the boids in `boid_state` according to their orientations and `paramSpeed`. If a boid moves past the edge of the box, it reappears on the opposite edge (periodic boundary conditions).

# Arguments
- `boid_state`: The current state of the boids.
- `paramSpeed`: The speed at which the boids move.
- `box_size`: A tuple representing the size of the box in each dimension.

# Example
```julia
movePeriodic!(boid_state, 1.0, (5, 5))
```
"""
function movePeriodic!(boid_state::BoidState, paramSpeed::Float64, box_size::Tuple)
    for i in eachindex(boid_state.positions)
        boid_state.positions[i] .+= boid_state.orientations[i] * paramSpeed
        for j in eachindex(boid_state.positions[i])
            if boid_state.positions[i][j] > box_size[j]
                boid_state.positions[i][j] -= box_size[j]
            end
            if boid_state.positions[i][j] < 0
                boid_state.positions[i][j] += box_size[j]
            end
        end
    end
end


function normReflect(v1::Vector{Float64}, v2::Vector{Float64})
    norm(v1 - v2)
end

"""
    normPeriodic(v1::Vector{Float64}, v2::Vector{Float64}, box_size::Tuple)

Calculate the distance between two "position" vectors `v1` and `v2` in a box with periodic boundary conditions.

# Arguments
- `v1`: The first vector.
- `v2`: The second vector.
- `box_size`: A tuple representing the size of the box in each dimension.

# Returns
- The minimum distance between `v1` and `v2` in the periodic box.

# Example
```julia
normPeriodic([1.0, 2.0], [3.0, 4.0], (5, 5))
```
"""
function normPeriodic(v1::Vector{Float64}, v2::Vector{Float64}, box_size::Tuple)
    temp = abs.(v1 - v2)
    temp = min.(temp, box_size .- temp)
    norm(temp)
end

"""
    runSimple2DBoids(numboids::Int64, steps::Int64, size::Float64, vel::Float64, paramMaxRad::Float64, paramMinRad::Float64, paramPosWeight::Float64, paramRepWeight::Float64, paramVelWeight::Float64)

Simulate the movement of `numboids` boids over `steps` time steps in a 2D box of size `size x size`.

The boids move with speed `vel` and interact with each other according to the given parameters.

# Arguments
- `numboids`: The number of boids.
- `steps`: The number of time steps.
- `size`: The size of the box in each dimension.
- `vel`: The speed at which the boids move.
- `paramMaxRad`: The maximum radius within which to consider other boids as neighbors.
- `paramMinRad`: The minimum radius within which other boids are considered too close and are repelled.
- `paramPosWeight`: The weight to give the effect of the average position of the neighbors when changing direction towards it.
- `paramRepWeight`: The weight to give to the repulsion from boids that are too close.
- `paramVelWeight`: The weight dictating how much a boid's orientation is affected by its neighbors' orientations.

# Returns
- A vector of `BoidState` objects representing the state of the boids at each time step.

# Example
```julia
runSimple2DBoids(10, 100, 5.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5)
```
"""
function runSimple2DBoids(numboids::Int64,steps::Int64,size::Float64,vel::Float64,paramMaxRad::Float64,paramMinRad::Float64,paramPosWeight::Float64,paramRepWeight::Float64,paramVelWeight::Float64)
    out=[BoidState(numboids,(size,size))]
    for i in 1:steps
        newstate = deepcopy(out[end])
        update!(newstate,state->ruleOrientations!(state,paramMaxRad,paramVelWeight,(v1,v2)->normPeriodic(v1,v2,(size,size))),state->rulePosition!(state,paramMaxRad,paramMinRad,paramPosWeight,paramRepWeight,(v1,v2)->normPeriodic(v1,v2,(size,size))),state->movePeriodic!(state,vel,(size,size)))
        push!(out,newstate)
    end
    return out
end

"""
    runSimple2DBoids(numboids::Int64, steps::Int64, size::Float64, vel::Float64, paramMaxRad::Float64, paramMinRad::Float64, paramPosWeight::Float64, paramRepWeight::Float64, paramVelWeight::Float64)

Simulate the movement of `numboids` boids over `steps` time steps in a 3D box of size `size x size x size`.

The boids move with speed `vel` and interact with each other according to the given parameters.

# Arguments
- `numboids`: The number of boids.
- `steps`: The number of time steps.
- `size`: The size of the box in each dimension.
- `vel`: The speed at which the boids move.
- `paramMaxRad`: The maximum radius within which to consider other boids as neighbors.
- `paramMinRad`: The minimum radius within which other boids are considered too close and are repelled.
- `paramPosWeight`: The weight to give the effect of the average position of the neighbors when changing direction towards it.
- `paramRepWeight`: The weight to give to the repulsion from boids that are too close.
- `paramVelWeight`: The weight dictating how much a boid's orientation is affected by its neighbors' orientations.

# Returns
- A vector of `BoidState` objects representing the state of the boids at each time step.

# Example
```julia
runSimple2DBoids(10, 100, 5.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5)
```
"""
function runSimple3DBoids(numboids::Int64,steps::Int64,size::Float64,vel::Float64,paramMaxRad::Float64,paramMinRad::Float64,paramPosWeight::Float64,paramRepWeight::Float64,paramVelWeight::Float64)
    out=[BoidState(numboids,(size,size,size))]
    for i in 1:steps
        newstate = deepcopy(out[end])
        update!(newstate,state->ruleOrientations!(state,paramMaxRad,paramVelWeight,(v1,v2)->normPeriodic(v1,v2,(size,size,size))),state->rulePosition!(state,paramMaxRad,paramMinRad,paramPosWeight,paramRepWeight,(v1,v2)->normPeriodic(v1,v2,(size,size,size))),state->movePeriodic!(state,vel,(size,size,size)))
        push!(out,newstate)
    end
    return out
end

end