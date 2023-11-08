# Boids

`Boids.jl` is a Julia package for simulating the behavior of boids, or bird-like objects, in a 2D space. It allows implementions of the classic boids algorithm(s), which models three basic behaviors: separation, alignment, and cohesion. 

## Installation

You can install `Boids.jl` using the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add https://github.com/jeanfdp/Boids.jl
```

## Usage
A simple implementation calculating the positions of boids in a periodic square box
```julia
using Boids

numboids=50
steps=200
size=1000.
vel=3.
paramMaxRad=100.
paramMinRad=10.
paramPosWeight=0.04
paramRepWeight=0.15
paramVelWeight=0.04

states=runSimple2DBoids(numboids,steps,size,vel,paramMaxRad,paramMinRad,paramPosWeight,paramRepWeight,paramVelWeight)
positions=[state.positions for state in states]
```



[![Build Status](https://github.com/jeanfdp/Boids.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/jeanfdp/Boids.jl/actions/workflows/CI.yml?query=branch%3Amaster)
