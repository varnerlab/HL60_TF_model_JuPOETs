include("Include.jl")

# Estimates model parameters for HL60 model -
function main()

  number_of_subdivisions = 10
  number_of_parameters = 153
  number_of_objectives = 12

  # Setup initial parameter guess -
  initial_parameter_array = zeros(number_of_parameters);

  # Load the data dictionary -
  data_dictionary = DataDictionary(0.0,0.0,0.0);

  # Update data dictionary to match new parameters before calculating obj
  parameter_mapping_array = data_dictionary["parameter_name_mapping_array"]
  for index = 1:length(parameter_mapping_array)
    if index <= data_dictionary["number_of_binding"]
      value = data_dictionary["binding_parameter_dictionary"][parameter_mapping_array[index]];
      initial_parameter_array[index] = value;
    else
      value = data_dictionary["control_parameter_dictionary"][parameter_mapping_array[index]];
      initial_parameter_array[index] = value;
    end
  end

  # Setup the error archives -
  ec_array = zeros(number_of_objectives)
  pc_array = zeros(number_of_parameters)
  for index in collect(1:number_of_subdivisions)

    # Run JuPOETs -
    (EC,PC,RA) = estimate_ensemble(objective_function,neighbor_function,acceptance_probability_function,cooling_function,initial_parameter_array;rank_cutoff=4,maximum_number_of_iterations=10,show_trace=true)

    # Package -
    ec_array = [ec_array EC]
    pc_array = [pc_array PC]

    # Take the *best* value from the current ensemble, refine it and go around again -
    total_error = sum(ec_array[:,2:end],1)

    # Which col is the min error?
    min_index = indmin(total_error)
    @show index,total_error[min_index]

    # dump to disk (in case we need to restart)
    # error cache -
    path_name = "./ensemble/ec_array.dat.$(index)"
    writedlm(path_name,ec_array)

    # parameter cache -
    path_name = "./ensemble/pc_array.dat.$(index)"
    writedlm(path_name,pc_array)

    # Refine the best solution -
    initial_parameter_array = pc_array[:,min_index].*(1+0.15*randn(number_of_parameters))
  end

  #return (ec_array,pc_array)
end

# Run -
main()
