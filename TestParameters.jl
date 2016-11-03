# Load -
include("Include.jl")

# Load the data dictionary -
data_dictionary = DataDictionary(0.0,0.0,0.0);

# Load the parameter array -
parameter_array = readdlm("./parameter_best_v1.dat")

# Update data dictionary to match new parameters before calculating obj
parameter_mapping_array = data_dictionary["parameter_name_mapping_array"]
for index = 1:length(parameter_mapping_array)
  if index <= data_dictionary["number_of_binding"]
    data_dictionary["binding_parameter_dictionary"][parameter_mapping_array[index]] = parameter_array[index]
  else
    data_dictionary["control_parameter_dictionary"][parameter_mapping_array[index]] = parameter_array[index]
  end
end

# Setup the time -
time_start = 0.0
time_stop = 120.0
time_step_size = 0.01

# Run the simulation -
(T,X) = Simulation(time_start,time_stop,time_step_size,data_dictionary);
