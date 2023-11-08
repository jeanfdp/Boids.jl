module Boids

using LinearAlgebra: normalize,normalize!,norm
using Distributions: MvNormal

export BoidState, BoidState, update!, ruleOrientations!, rulePosition!, moveReflect!, movePeriodic!, normReflect, normPeriodic, runSimple2DBoids, runSimple3DBoids

struct BoidState
    positions::Vector{Vector{Float64}}
    orientations::Vector{Vector{Float64}}
end

function BoidState(n::Int, box_size::Tuple)
    dim = length(box_size)
    positions = [rand(dim) .* box_size for _ in 1:n]
    orientations = [normalize(rand(MvNormal(dim, 1))) for _ in 1:n]
    BoidState(positions, orientations)
end

function update!(state::BoidState,ruleOrientations!::Function,rulePositions!::Function,move!::Function)
    ruleOrientations!(state)
    rulePositions!(state)
    move!(state)
end


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
                temprepulse .+= boid_state.positions[j]
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

function normPeriodic(v1::Vector{Float64}, v2::Vector{Float64}, box_size::Tuple)
    temp = abs.(v1 - v2)
    temp = min.(temp, box_size .- temp)
    norm(temp)
end

function runSimple2DBoids(numboids::Int64,steps::Int64,size::Float64,vel::Float64,paramMaxRad::Float64,paramMinRad::Float64,paramPosWeight::Float64,paramRepWeight::Float64,paramVelWeight::Float64)
    out=[BoidState(numboids,(size,size))]
    for i in 1:steps
        newstate = deepcopy(out[end])
        update!(newstate,state->ruleOrientations!(state,paramMaxRad,paramVelWeight,(v1,v2)->normPeriodic(v1,v2,(size,size))),state->rulePosition!(state,paramMaxRad,paramMinRad,paramPosWeight,paramRepWeight,(v1,v2)->normPeriodic(v1,v2,(size,size))),state->movePeriodic!(state,vel,(size,size)))
        push!(out,newstate)
    end
    return out
end

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