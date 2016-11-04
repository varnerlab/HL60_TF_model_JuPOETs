# include -
include("Include.jl")

# How many parameters?
number_of_parameters = 153;
number_of_initial_sets = 5;
number_of_final_sets = 200;

# Simple script to generate ensemble by sampling -
tmp_parameter_array = zeros(number_of_parameters,number_of_initial_sets);
for parameter_set_index = 1:number_of_initial_sets

  # Setup tmp filename -
  tmp_filename = "./parameter_best_v1.dat."*string(parameter_set_index);
  parameter_set = readdlm(tmp_filename);

  for parameter_index = 1:number_of_parameters
    tmp_parameter_array[parameter_index,parameter_set_index] = parameter_set[parameter_index];
  end
end

# Establish bounds -
bounds_array = zeros(number_of_parameters,2);
for parameter_index = 1:number_of_parameters

  # min and max -
  minimum_value = minimum(tmp_parameter_array[parameter_index,:]);
  maximum_value = maximum(tmp_parameter_array[parameter_index,:]);

  # perturbation -
  perturbation_size = 1;
  bounds_array[parameter_index,1] = (1/perturbation_size)*minimum_value;
  bounds_array[parameter_index,2] = perturbation_size*maximum_value;
end

# Generate an ensemble of parameter sets -
ensemble_array = zeros(number_of_parameters,number_of_final_sets);
sobol_sequence = SobolSeq(number_of_parameters,bounds_array[:,1],bounds_array[:,2]);
for parameter_set_index = 1:number_of_final_sets

  # sample the parameter array -
  parameter_array = next(sobol_sequence)

  for parameter_index = 1:number_of_parameters
    ensemble_array[parameter_index,parameter_set_index] = parameter_array[parameter_index];
  end
end

# Solve the model to calculate the error -
error_array = zeros(number_of_final_sets);
for parameter_set_index = 1:number_of_final_sets

  # sample the ensemble_array -
  parameter_array = ensemble_array[:,parameter_set_index]

  # error -
  error_value = sum(objective_function(parameter_array));
  error_array[parameter_set_index] = error_value;

  @show (parameter_set_index,error_value)
end

# Filter -
minimum_error = minimum(error_array);
variance_error = var(error_array);
epsilon = 1.0;
error_cutoff = minimum_error*(1.0+epsilon);
ensemble_index_array = Int[];
for parameter_set_index = 1:number_of_final_sets

  error_value = error_array[parameter_set_index];

  if (error_value<=error_cutoff)
    push!(ensemble_index_array,parameter_set_index);
  else

    # can we accept this set?
    error_difference = (error_value - minimum_error)^2
    prob_keep_set = exp(-1*error_difference/variance_error);
    check_role = rand();
    if (prob_keep_set>check_role)
      push!(ensemble_index_array,parameter_set_index);
    end
  end
end

# Last ... run the model w/these sets, save the data -
for (local_index,parameter_set_index) in enumerate(ensemble_index_array)

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
  data_filename = "./ensemble/sim_data_EI_"*string(parameter_set_index)*"_LI_"*string(local_index)*".dat"
  writedlm(data_filename,local_data);
end

# Write the simple ensemble -
filename = "./simple_ensemble_parameter_sets.dat"
writedlm(filename,ensemble_array);
