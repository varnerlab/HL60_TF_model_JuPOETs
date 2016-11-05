# include -
include("Include.jl")

# setup the calculation -
number_of_samples = 5
number_of_parameters = 153

# load up the parameters sets -
ensemble_array = zeros(number_of_parameters,1);
for sample_index = 1:number_of_samples

  filename = "./parameter_best_v1.dat."*string(sample_index)
  local_array = readdlm(filename)

  ensemble_array = hcat(ensemble_array,local_array);
end

# cut leading col -
ensemble_array = ensemble_array[:,2:end]

# Last ... run the model w/these sets, save the data -
for parameter_set_index = 1:number_of_samples

  # sample the ensemble_array -
  parameter_array = ensemble_array[:,parameter_set_index]

  # Load the data dictionary -
  data_dictionary = DataDictionary(0.0,0.0,0.0);

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

  # dump data to disk -
  local_data = [T X];
  data_filename = "./best_ensemble/sim_data_PSI_"*string(parameter_set_index)*".dat"
  writedlm(data_filename,local_data);
end
