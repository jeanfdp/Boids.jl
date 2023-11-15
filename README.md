# Boids

`Boids.jl` is a Julia package for simulating the behavior of boids, or bird-like objects, in both 2D and 3D spaces. This package allows you to implement and explore the classic boids algorithm(s), which model three fundamental behaviors: separation, alignment, and cohesion. It is intended to allow a collection of different implementation of separation, alignment and cohesion rules (as well as potentially more complicated and exotic rules as well). Contributions are more than welcome and discussed in the 'Contribution' section below.


## Installation

You can easily install `Boids.jl` using the Julia package manager. Open the Julia REPL, and type the following commands:

```julia
using Pkg
Pkg.add("Boids")
```


## Basic Usage

### Creating a BoidState
To work with the Boids package, you first need to create a BoidState object. This object represents the state of a flock of boids in some space. You can either initialize it with specific positions and orientations or create random initial states within a specified box.

```julia
using Boids

# Initialize a BoidState with specific positions and orientations
state = BoidState([[1.0, 2.0], [3.0, 4.0]], [[0.5, 0.5], [0.7, 0.3]])

# Create a random BoidState with 'n' boids in a 2D box of size 'box_size'
n = 10
box_size = (5.0, 5.0)
state = BoidState(n, box_size)

# Create a random BoidState with 'n' boids in a 3D box of size 'box_size'
n = 10
box_size = (5.0, 5.0,5.0)
state = BoidState(n, box_size)
```

### Simulating Boid Behaviour
In the most general case you can simulate the movement and interaction of boids over a series of time steps using the update! function. This function applies rules for orientations and positions and moves the boids accordingly. You can specify your custom rules and movement functions.
```julia
# Define custom rules and movement functions
function customRuleOrientations!(boid_state::BoidState)
    # Your custom orientation rule implementation here
end

function customRulePositions!(boid_state::BoidState)
    # Your custom position rule implementation here
end

function customMove!(boid_state::BoidState)
    # Your custom movement implementation here
end

# Apply the custom rules and movement
update!(state, customRuleOrientations!, customRulePositions!, customMove!)
```
This package is intended to eventually include whatever rules people find useful, and to gather them in a coherent package.

### Example Simulation
We can use some of the built in update rules in order to run a simple 2D simulation in a periodic square box:
```julia
numboids = 50
steps = 200
size = 1000.0
vel = 3.0
paramMaxRad = 100.0
paramMinRad = 10.0
paramPosWeight = 0.04
paramRepWeight = 0.15
paramVelWeight = 0.04

# Simulate the boid behavior
states = runSimple2DBoids(numboids, steps, size, vel, paramMaxRad, paramMinRad, paramPosWeight, paramRepWeight, paramVelWeight)

# Extract positions from the simulation states
positions = [state.positions for state in states]
```


## Contributing

If you're interested in contributing to `Boids.jl`, we welcome your contributions and appreciate your help in making the package better. To get started, follow these guidelines:

1. **Fork the Repository:** Go to the [GitHub repository](https://github.com/jeanfdp/Boids.jl) and click the "Fork" button to create your own fork of the project.

2. **Clone the Repository:** Clone your fork of the repository to your local machine.

3. **Create a New Branch:** Create a new branch with a descriptive name for your feature, enhancement, or bug fix.

4. **Make Changes:** Make your changes and ensure they're well-tested.

5. **Commit Changes:** Commit your changes with clear and concise commit messages.

6. **Push Changes:** Push your branch to your GitHub fork.

7. **Submit a Pull Request:** Create a pull request from your branch to the main repository, describing the purpose of your PR and the changes you've made.

8. **Participate in Review:** Participate in the review process and address any feedback if required.

We appreciate your contributions and look forward to collaborating with you to improve `Boids.jl`.


## License

`Boids.jl` is open-source and available under the [MIT License](https://github.com/jeanfdp/Boids.jl/blob/master/LICENSE).




[![Build Status](https://github.com/jeanfdp/Boids.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/jeanfdp/Boids.jl/actions/workflows/CI.yml?query=branch%3Amaster)
