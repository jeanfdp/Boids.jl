module Boids

using LinearAlgebra: normalize
using Distributions: MvNormal

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

function update!(state::BoidState, box_size::Tuple, v::Float64, w::Float64,ruleOrientations!::Function,rulePositions!::Function)
    dim = length(box_size)
    n = length(state.positions)
    for i in 1:n
        ruleOrientations!(state, i, w)
        rulePositions!(state, i, v, box_size)
    end
end

end
