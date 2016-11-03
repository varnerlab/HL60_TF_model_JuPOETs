# includes -
include("Kinetics.jl")
include("Control.jl")
include("Inputs.jl")
include("DataDictionary.jl")
include("SolveBalances.jl")
include("Balances.jl")
using ODE

# Required by the tf_network_lib -
using JSON
using POETs
include("Error.jl")
include("Simulation.jl")
include("tf_network_lib.jl")

# Required by the error function -
using PyCall
@pyimport numpy as np

# Plotting -
using PyPlot
