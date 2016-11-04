function process_ensemble_data(time_array,species_index)

  # Load the data from disk -
  path_to_sim_files = "./ensemble"
  file_extension = ".dat"
  searchdir(path,key) = filter(x->contains(x,key),readdir(path))

  # build src file list -
  list_of_sim_files = searchdir(path_to_sim_files,file_extension);

  # initialize data array -
  raw_data_array = zeros(length(time_array),1);

  # go thru the src file list, and copy the files to the output path -
  for (sim_file_index,sim_file) in enumerate(list_of_sim_files)

    # path to sim file -
    sim_file_path = "./ensemble/"*string(sim_file)

    # load the data array -
    simulation_data_array = readdlm(sim_file_path);

    # Get the time -
    simulation_time_array = simulation_data_array[:,1];
    simulation_state_array = simulation_data_array[:,species_index+1];

    # interpolate -
    interpolated_data_array = np.interp(time_array,simulation_time_array,simulation_state_array);

    @show (sim_file_index,sim_file)
    raw_data_array = hcat(raw_data_array,transpose(interpolated_data_array))
  end

  # Generate the scaled data -
  raw_data_array = raw_data_array[:,2:end]
  (number_of_rows,number_of_cols) = size(raw_data_array);
  scaled_data_array = zeros(number_of_rows,number_of_cols);
  for col_index = 1:number_of_cols

    # grab the first elements -
    scale_element = raw_data_array[1,col_index];
    for row_index = 1:number_of_rows
      scaled_value_element = raw_data_array[row_index,col_index]*(1/scale_element);
      scaled_data_array[row_index,col_index] = scaled_value_element;
    end
  end

  return (raw_data_array,scaled_data_array)
end
