# include -
include("Include.jl")

# What index do I want to see?
plot_index = 41;

# Load the simulation result set -
result_set = process_simulation_files("./simple_ensemble")
for (index,data_array) in enumerate(result_set)

  @show index

  time_array = data_array[:,1];
  state_array = data_array[:,plot_index+1];

  plot(time_array,state_array);
end
