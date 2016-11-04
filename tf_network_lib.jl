# ================ INITIALIZATION CODE DO NOT EDIT ====================================== %
# some global parameters -
BIG = 1e10
SMALL = 1e-8

# Load the experiment specifications -
tmp_value = JSON.parsefile("./experiments/Experiments.json")
experiment_array = tmp_value["experiment_array"]

# preload the data -
cached_data_dictionary = Dict()
for (experiment_index,experiment_dictionary) in enumerate(experiment_array)

  # Grab the data file -
  data_file_path = experiment_dictionary["data_file"]
  experiment_id = experiment_dictionary["experiment_id"]

  # Load the data -
  experimental_data_array = readdlm(data_file_path)

  # Cache the data w/experiment_id -
  cached_data_dictionary[experiment_id] = experimental_data_array
end

# preload the error function calls -
cached_error_function_dictionary = Dict()
for (experiment_index,experiment_dictionary) in enumerate(experiment_array)

  # Grab the experimental id -
  experiment_id = experiment_dictionary["experiment_id"]

  # Grab the error functions from the experiment_dictionary -
  error_function = eval(parse(experiment_dictionary["error_function"]))

  # cache these functions -
  cached_error_function_dictionary[experiment_id] = error_function
end
# ================ INITIALIZATION CODE DO NOT EDIT ====================================== %


function objective_function(parameter_array)

  # Script to solve the balance equations -
  time_start = 0.0
  time_stop = 120.0
  time_step_size = 0.01

  # Load the data dictionary -
  data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

  # Update data dictionary to match new parameters before calculating obj
  parameter_mapping_array = data_dictionary["parameter_name_mapping_array"]
  for index = 1:length(parameter_mapping_array)
    if index <= data_dictionary["number_of_binding"]
      data_dictionary["binding_parameter_dictionary"][parameter_mapping_array[index]] = parameter_array[index]
    else
      data_dictionary["control_parameter_dictionary"][parameter_mapping_array[index]] = parameter_array[index]
    end
  end

  # how many objectives do we have?
  obj_array = BIG*ones(13,1)

  # Call simulation routine -
  # (run the model to SS, and then set the ICs to the SS for this parameter set)
  (time_array,simulation_state_array) = Simulation(time_start,time_stop,time_step_size,data_dictionary);

  # Call the error functions -
  # loop through the experimental dictionary, and call the appropriate error function -
  for (experiment_index,experiment_dictionary) in enumerate(experiment_array)

    # Get the error function pointer -
    error_function_pointer = experiment_dictionary["error_function"]
    output_index = parse(Int,experiment_dictionary["output_index"])
    experiment_id = experiment_dictionary["experiment_id"]
    species_symbol = experiment_dictionary["protein_symbol"];

    # Get the experimental data array -
    experimental_data_array = cached_data_dictionary[experiment_id]

    # Call the error function -
    error_function = cached_error_function_dictionary[experiment_id]
    error_value = error_function(experimental_data_array,time_array,simulation_state_array,output_index,species_symbol,data_dictionary)

    # Add the error to the objective array -
    obj_array[experiment_index] = error_value;
  end

  # @show obj_array

  # return the error array -
  return obj_array;
end

function neighbor_function(parameter_array)

  SIGMA = 0.01
  number_of_parameters = length(parameter_array)

  # calculate new parameters -
  perturbed_parameter_array = parameter_array.*(1+SIGMA*randn(number_of_parameters));

  # Parameters -
  # "n_gene_AP1_gene_AhR"	;	# 1
  # "K_gene_AP1_gene_AhR"	;	# 2
  # "n_gene_AP1_gene_PU1"	;	# 3
  # "K_gene_AP1_gene_PU1"	;	# 4
  # "n_gene_AP1_gene_PPARg"	;	# 5
  # "K_gene_AP1_gene_PPARg"	;	# 6
  # "n_gene_AhR_gene_Trigger"	;	# 7
  # "K_gene_AhR_gene_Trigger"	;	# 8
  # "n_gene_CD11b_gene_PU1_gene_cRAF"	;	# 9
  # "K_gene_CD11b_gene_PU1_gene_cRAF"	;	# 10
  # "n_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF"	;	# 11
  # "K_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF"	;	# 12
  # "n_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF"	;	# 13
  # "K_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF"	;	# 14
  # "n_gene_CEBPa_gene_Trigger"	;	# 15
  # "K_gene_CEBPa_gene_Trigger"	;	# 16
  # "n_gene_CEBPa_gene_PPARg"	;	# 17
  # "K_gene_CEBPa_gene_PPARg"	;	# 18
  # "n_gene_CEBPa_gene_CEBPa"	;	# 19
  # "K_gene_CEBPa_gene_CEBPa"	;	# 20
  # "n_gene_CEBPa_gene_GFI1"	;	# 21
  # "K_gene_CEBPa_gene_GFI1"	;	# 22
  # "n_gene_E2F_gene_E2F"	;	# 23
  # "K_gene_E2F_gene_E2F"	;	# 24
  # "n_gene_E2F_gene_PPARg"	;	# 25
  # "K_gene_E2F_gene_PPARg"	;	# 26
  # "n_gene_E2F_gene_CEBPa"	;	# 27
  # "K_gene_E2F_gene_CEBPa"	;	# 28
  # "n_gene_E2F_gene_GFI1"	;	# 29
  # "K_gene_E2F_gene_GFI1"	;	# 30
  # "n_gene_E2F_gene_cRAF"	;	# 31
  # "K_gene_E2F_gene_cRAF"	;	# 32
  # "n_gene_EGR1_gene_Trigger"	;	# 33
  # "K_gene_EGR1_gene_Trigger"	;	# 34
  # "n_gene_EGR1_gene_PU1"	;	# 35
  # "K_gene_EGR1_gene_PU1"	;	# 36
  # "n_gene_EGR1_gene_PPARg"	;	# 37
  # "K_gene_EGR1_gene_PPARg"	;	# 38
  # "n_gene_EGR1_gene_GFI1"	;	# 39
  # "K_gene_EGR1_gene_GFI1"	;	# 40
  # "n_gene_GFI1_gene_CEBPa"	;	# 41
  # "K_gene_GFI1_gene_CEBPa"	;	# 42
  # "n_gene_GFI1_gene_EGR1"	;	# 43
  # "K_gene_GFI1_gene_EGR1"	;	# 44
  # "n_gene_IRF1_gene_Trigger"	;	# 45
  # "K_gene_IRF1_gene_Trigger"	;	# 46
  # "n_gene_IRF1_gene_AhR"	;	# 47
  # "K_gene_IRF1_gene_AhR"	;	# 48
  # "n_gene_IRF1_gene_PPARg"	;	# 49
  # "K_gene_IRF1_gene_PPARg"	;	# 50
  # "n_gene_OCT1_gene_PPARg"	;	# 51
  # "K_gene_OCT1_gene_PPARg"	;	# 52
  # "n_gene_OCT4_gene_Trigger"	;	# 53
  # "K_gene_OCT4_gene_Trigger"	;	# 54
  # "n_gene_OCT4_gene_AhR"	;	# 55
  # "K_gene_OCT4_gene_AhR"	;	# 56
  # "n_gene_OCT4_gene_cRAF"	;	# 57
  # "K_gene_OCT4_gene_cRAF"	;	# 58
  # "n_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF"	;	# 59
  # "K_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF"	;	# 60
  # "n_gene_P21_gene_GFI1"	;	# 61
  # "K_gene_P21_gene_GFI1"	;	# 62
  # "n_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF"	;	# 63
  # "K_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF"	;	# 64
  # "n_gene_P47Phox_gene_PPARg"	;	# 65
  # "K_gene_P47Phox_gene_PPARg"	;	# 66
  # "n_gene_PPARg_gene_Trigger"	;	# 67
  # "K_gene_PPARg_gene_Trigger"	;	# 68
  # "n_gene_PPARg_gene_CEBPa"	;	# 69
  # "K_gene_PPARg_gene_CEBPa"	;	# 70
  # "n_gene_PPARg_gene_EGR1"	;	# 71
  # "K_gene_PPARg_gene_EGR1"	;	# 72
  # "n_gene_PPARg_gene_PU1"	;	# 73
  # "K_gene_PPARg_gene_PU1"	;	# 74
  # "n_gene_PPARg_gene_AP1"	;	# 75
  # "K_gene_PPARg_gene_AP1"	;	# 76
  # "n_gene_PU1_gene_Trigger"	;	# 77
  # "K_gene_PU1_gene_Trigger"	;	# 78
  # "n_gene_PU1_gene_CEBPa"	;	# 79
  # "K_gene_PU1_gene_CEBPa"	;	# 80
  # "n_gene_PU1_gene_PU1"	;	# 81
  # "K_gene_PU1_gene_PU1"	;	# 82
  # "n_gene_PU1_gene_AP1"	;	# 83
  # "K_gene_PU1_gene_AP1"	;	# 84
  # "n_gene_PU1_gene_OCT1"	;	# 85
  # "K_gene_PU1_gene_OCT1"	;	# 86
  # "n_gene_PU1_gene_AhR"	;	# 87
  # "K_gene_PU1_gene_AhR"	;	# 88
  # "n_gene_PU1_gene_GFI1"	;	# 89
  # "K_gene_PU1_gene_GFI1"	;	# 90
  # "W_gene_AP1_RNAP"	;	# 91
  # "W_gene_AP1_gene_AhR"	;	# 92
  # "W_gene_AP1_gene_PU1"	;	# 93
  # "W_gene_AP1_gene_PPARg"	;	# 94
  # "W_gene_AhR_RNAP"	;	# 95
  # "W_gene_AhR_gene_Trigger"	;	# 96
  # "W_gene_CD11b_RNAP"	;	# 97
  # "W_gene_CD11b_gene_PU1_gene_cRAF"	;	# 98
  # "W_gene_CD14_RNAP"	;	# 99
  # "W_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF"	;	# 100
  # "W_gene_CD38_RNAP"	;	# 101
  # "W_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF"	;	# 102
  # "W_gene_CEBPa_RNAP"	;	# 103
  # "W_gene_CEBPa_gene_Trigger"	;	# 104
  # "W_gene_CEBPa_gene_PPARg"	;	# 105
  # "W_gene_CEBPa_gene_CEBPa"	;	# 106
  # "W_gene_CEBPa_gene_GFI1"	;	# 107
  # "W_gene_E2F_RNAP"	;	# 108
  # "W_gene_E2F_gene_E2F"	;	# 109
  # "W_gene_E2F_gene_PPARg"	;	# 110
  # "W_gene_E2F_gene_CEBPa"	;	# 111
  # "W_gene_E2F_gene_GFI1"	;	# 112
  # "W_gene_E2F_gene_cRAF"	;	# 113
  # "W_gene_EGR1_RNAP"	;	# 114
  # "W_gene_EGR1_gene_Trigger"	;	# 115
  # "W_gene_EGR1_gene_PU1"	;	# 116
  # "W_gene_EGR1_gene_PPARg"	;	# 117
  # "W_gene_EGR1_gene_GFI1"	;	# 118
  # "W_gene_GFI1_RNAP"	;	# 119
  # "W_gene_GFI1_gene_CEBPa"	;	# 120
  # "W_gene_GFI1_gene_EGR1"	;	# 121
  # "W_gene_IRF1_RNAP"	;	# 122
  # "W_gene_IRF1_gene_Trigger"	;	# 123
  # "W_gene_IRF1_gene_AhR"	;	# 124
  # "W_gene_IRF1_gene_PPARg"	;	# 125
  # "W_gene_OCT1_RNAP"	;	# 126
  # "W_gene_OCT1_gene_PPARg"	;	# 127
  # "W_gene_OCT4_RNAP"	;	# 128
  # "W_gene_OCT4_gene_Trigger"	;	# 129
  # "W_gene_OCT4_gene_AhR"	;	# 130
  # "W_gene_OCT4_gene_cRAF"	;	# 131
  # "W_gene_P21_RNAP"	;	# 132
  # "W_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF"	;	# 133
  # "W_gene_P21_gene_GFI1"	;	# 134
  # "W_gene_P47Phox_RNAP"	;	# 135
  # "W_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF"	;	# 136
  # "W_gene_P47Phox_gene_PPARg"	;	# 137
  # "W_gene_PPARg_RNAP"	;	# 138
  # "W_gene_PPARg_gene_Trigger"	;	# 139
  # "W_gene_PPARg_gene_CEBPa"	;	# 140
  # "W_gene_PPARg_gene_EGR1"	;	# 141
  # "W_gene_PPARg_gene_PU1"	;	# 142
  # "W_gene_PPARg_gene_AP1"	;	# 143
  # "W_gene_PU1_RNAP"	;	# 144
  # "W_gene_PU1_gene_Trigger"	;	# 145
  # "W_gene_PU1_gene_CEBPa"	;	# 146
  # "W_gene_PU1_gene_PU1"	;	# 147
  # "W_gene_PU1_gene_AP1"	;	# 148
  # "W_gene_PU1_gene_OCT1"	;	# 149
  # "W_gene_PU1_gene_AhR"	;	# 150
  # "W_gene_PU1_gene_GFI1"	;	# 151
  # "W_gene_Trigger_RNAP"	;	# 152
  # "W_gene_cRAF_RNAP"	;	# 153

  # Setup my upper bound, and lower bounds on parameters -
  lower_bound_array = zeros(number_of_parameters)
  upper_bound_array = [
    4.0   ;   # "n_gene_AP1_gene_AhR"	;	# 1
    500.0 ;   # "K_gene_AP1_gene_AhR"	;	# 2
    4.0   ;   # "n_gene_AP1_gene_PU1"	;	# 3
    500.0 ;   # "K_gene_AP1_gene_PU1"	;	# 4
    4.0   ;   # "n_gene_AP1_gene_PPARg"	;	# 5
    500.0 ;   # "K_gene_AP1_gene_PPARg"	;	# 6
    4.0   ;   # "n_gene_AhR_gene_Trigger"	;	# 7
    500.0 ;   # "K_gene_AhR_gene_Trigger"	;	# 8
    4.0   ;  # "n_gene_CD11b_gene_PU1_gene_cRAF"	;	# 9
    500.0 ;  # "K_gene_CD11b_gene_PU1_gene_cRAF"	;	# 10
    4.0   ;  # "n_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF"	;	# 11
    500.0 ;  # "K_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF"	;	# 12
    4.0   ;  # "n_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF"	;	# 13
    500.0 ;  # "K_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF"	;	# 14
    4.0   ;  # "n_gene_CEBPa_gene_Trigger"	;	# 15
    500.0 ;  # "K_gene_CEBPa_gene_Trigger"	;	# 16
    4.0   ;  # "n_gene_CEBPa_gene_PPARg"	;	# 17
    500.0 ;  # "K_gene_CEBPa_gene_PPARg"	;	# 18
    4.0   ;  # "n_gene_CEBPa_gene_CEBPa"	;	# 19
    500.0 ;  # "K_gene_CEBPa_gene_CEBPa"	;	# 20
    4.0   ;  # "n_gene_CEBPa_gene_GFI1"	;	# 21
    500.0 ;  # "K_gene_CEBPa_gene_GFI1"	;	# 22
    4.0   ;  # "n_gene_E2F_gene_E2F"	;	# 23
    500.0 ;  # "K_gene_E2F_gene_E2F"	;	# 24
    4.0   ;  # "n_gene_E2F_gene_PPARg"	;	# 25
    500.0 ;  # "K_gene_E2F_gene_PPARg"	;	# 26
    4.0   ;  # "n_gene_E2F_gene_CEBPa"	;	# 27
    500.0 ;  # "K_gene_E2F_gene_CEBPa"	;	# 28
    4.0   ;  # "n_gene_E2F_gene_GFI1"	;	# 29
    500.0 ;  # "K_gene_E2F_gene_GFI1"	;	# 30
    4.0   ;  # "n_gene_E2F_gene_cRAF"	;	# 31
    500.0 ;  # "K_gene_E2F_gene_cRAF"	;	# 32
    4.0   ;  # "n_gene_EGR1_gene_Trigger"	;	# 33
    500.0 ;  # "K_gene_EGR1_gene_Trigger"	;	# 34
    4.0   ;  # "n_gene_EGR1_gene_PU1"	;	# 35
    500.0 ;  # "K_gene_EGR1_gene_PU1"	;	# 36
    4.0   ;  # "n_gene_EGR1_gene_PPARg"	;	# 37
    500.0 ;  # "K_gene_EGR1_gene_PPARg"	;	# 38
    4.0   ;  # "n_gene_EGR1_gene_GFI1"	;	# 39
    500.0 ;  # "K_gene_EGR1_gene_GFI1"	;	# 40
    4.0   ;  # "n_gene_GFI1_gene_CEBPa"	;	# 41
    500.0 ;  # "K_gene_GFI1_gene_CEBPa"	;	# 42
    4.0   ;  # "n_gene_GFI1_gene_EGR1"	;	# 43
    500.0 ;  # "K_gene_GFI1_gene_EGR1"	;	# 44
    4.0   ;  # "n_gene_IRF1_gene_Trigger"	;	# 45
    500.0 ;  # "K_gene_IRF1_gene_Trigger"	;	# 46
    4.0   ;  # "n_gene_IRF1_gene_AhR"	;	# 47
    500.0 ;  # "K_gene_IRF1_gene_AhR"	;	# 48
    4.0   ;  # "n_gene_IRF1_gene_PPARg"	;	# 49
    500.0 ;  # "K_gene_IRF1_gene_PPARg"	;	# 50
    4.0   ;  # "n_gene_OCT1_gene_PPARg"	;	# 51
    500.0 ;  # "K_gene_OCT1_gene_PPARg"	;	# 52
    4.0   ;  # "n_gene_OCT4_gene_Trigger"	;	# 53
    500.0 ;  # "K_gene_OCT4_gene_Trigger"	;	# 54
    4.0   ;  # "n_gene_OCT4_gene_AhR"	;	# 55
    500.0 ;  # "K_gene_OCT4_gene_AhR"	;	# 56
    4.0   ;  # "n_gene_OCT4_gene_cRAF"	;	# 57
    500.0 ;  # "K_gene_OCT4_gene_cRAF"	;	# 58
    4.0   ;  # "n_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF"	;	# 59
    500.0 ;  # "K_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF"	;	# 60
    4.0   ;  # "n_gene_P21_gene_GFI1"	;	# 61
    500.0 ;  # "K_gene_P21_gene_GFI1"	;	# 62
    4.0   ;  # "n_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF"	;	# 63
    500.0 ;  # "K_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF"	;	# 64
    4.0   ;  # "n_gene_P47Phox_gene_PPARg"	;	# 65
    500.0 ;  # "K_gene_P47Phox_gene_PPARg"	;	# 66
    4.0   ;  # "n_gene_PPARg_gene_Trigger"	;	# 67
    500.0 ;  # "K_gene_PPARg_gene_Trigger"	;	# 68
    4.0   ;  # "n_gene_PPARg_gene_CEBPa"	;	# 69
    500.0 ;  # "K_gene_PPARg_gene_CEBPa"	;	# 70
    4.0   ;  # "n_gene_PPARg_gene_EGR1"	;	# 71
    500.0 ;  # "K_gene_PPARg_gene_EGR1"	;	# 72
    4.0   ;  # "n_gene_PPARg_gene_PU1"	;	# 73
    500.0 ;  # "K_gene_PPARg_gene_PU1"	;	# 74
    4.0   ;  # "n_gene_PPARg_gene_AP1"	;	# 75
    500.0 ;  # "K_gene_PPARg_gene_AP1"	;	# 76
    4.0   ;  # "n_gene_PU1_gene_Trigger"	;	# 77
    500.0 ;  # "K_gene_PU1_gene_Trigger"	;	# 78
    4.0   ;  # "n_gene_PU1_gene_CEBPa"	;	# 79
    500.0 ;  # "K_gene_PU1_gene_CEBPa"	;	# 80
    4.0   ;  # "n_gene_PU1_gene_PU1"	;	# 81
    500.0 ;  # "K_gene_PU1_gene_PU1"	;	# 82
    4.0   ;  # "n_gene_PU1_gene_AP1"	;	# 83
    500.0 ;  # "K_gene_PU1_gene_AP1"	;	# 84
    4.0   ;  # "n_gene_PU1_gene_OCT1"	;	# 85
    500.0 ;  # "K_gene_PU1_gene_OCT1"	;	# 86
    4.0   ;  # "n_gene_PU1_gene_AhR"	;	# 87
    500.0 ;  # "K_gene_PU1_gene_AhR"	;	# 88
    4.0   ;  # "n_gene_PU1_gene_GFI1"	;	# 89
    500.0 ;  # "K_gene_PU1_gene_GFI1"	;	# 90
    100.0 ; # "W_gene_AP1_RNAP"	;	# 91
    100.0 ;  # "W_gene_AP1_gene_AhR"	;	# 92
    100.0 ;  # "W_gene_AP1_gene_PU1"	;	# 93
    100.0 ;  # "W_gene_AP1_gene_PPARg"	;	# 94
    100.0 ;  # "W_gene_AhR_RNAP"	;	# 95
    100.0 ;  # "W_gene_AhR_gene_Trigger"	;	# 96
    100.0 ;  # "W_gene_CD11b_RNAP"	;	# 97
    100.0 ;  # "W_gene_CD11b_gene_PU1_gene_cRAF"	;	# 98
    100.0 ;  # "W_gene_CD14_RNAP"	;	# 99
    0.0   ;  # "W_gene_CD14_gene_PPARg_gene_CEBPa_gene_EGR1_gene_cRAF"	;	# 100
    100.0 ;  # "W_gene_CD38_RNAP"	;	# 101
    100.0 ;  # "W_gene_CD38_gene_IRF1_gene_PPARg_gene_Trigger_gene_cRAF"	;	# 102
    100.0 ;  # "W_gene_CEBPa_RNAP"	;	# 103
    100.0 ;  # "W_gene_CEBPa_gene_Trigger"	;	# 104
    100.0 ;  # "W_gene_CEBPa_gene_PPARg"	;	# 105
    100.0 ;  # "W_gene_CEBPa_gene_CEBPa"	;	# 106
    100.0 ;  # "W_gene_CEBPa_gene_GFI1"	;	# 107
    100.0 ;  # "W_gene_E2F_RNAP"	;	# 108
    100.0 ;  # "W_gene_E2F_gene_E2F"	;	# 109
    100.0 ;  # "W_gene_E2F_gene_PPARg"	;	# 110
    100.0 ;  # "W_gene_E2F_gene_CEBPa"	;	# 111
    100.0 ;  # "W_gene_E2F_gene_GFI1"	;	# 112
    100.0 ;  # "W_gene_E2F_gene_cRAF"	;	# 113
    100.0 ;  # "W_gene_EGR1_RNAP"	;	# 114
    100.0 ;  # "W_gene_EGR1_gene_Trigger"	;	# 115
    100.0 ;  # "W_gene_EGR1_gene_PU1"	;	# 116
    100.0 ;  # "W_gene_EGR1_gene_PPARg"	;	# 117
    100.0 ;  # "W_gene_EGR1_gene_GFI1"	;	# 118
    100.0 ;  # "W_gene_GFI1_RNAP"	;	# 119
    100.0 ;  # "W_gene_GFI1_gene_CEBPa"	;	# 120
    100.0 ;  # "W_gene_GFI1_gene_EGR1"	;	# 121
    100.0 ;  # "W_gene_IRF1_RNAP"	;	# 122
    100.0 ;  # "W_gene_IRF1_gene_Trigger"	;	# 123
    100.0 ;  # "W_gene_IRF1_gene_AhR"	;	# 124
    100.0 ;  # "W_gene_IRF1_gene_PPARg"	;	# 125
    100.0 ;  # "W_gene_OCT1_RNAP"	;	# 126
    100.0 ;  # "W_gene_OCT1_gene_PPARg"	;	# 127
    100.0 ;  # "W_gene_OCT4_RNAP"	;	# 128
    100.0 ;  # "W_gene_OCT4_gene_Trigger"	;	# 129
    100.0 ;  # "W_gene_OCT4_gene_AhR"	;	# 130
    100.0 ;  # "W_gene_OCT4_gene_cRAF"	;	# 131
    100.0 ;  # "W_gene_P21_RNAP"	;	# 132
    100.0 ;  # "W_gene_P21_gene_Trigger_gene_AP1_gene_PPARg_gene_PU1_gene_IRF1_gene_CEBPa_gene_cRAF"	;	# 133
    100.0 ;  # "W_gene_P21_gene_GFI1"	;	# 134
    100.0 ;  # "W_gene_P47Phox_RNAP"	;	# 135
    100.0 ;  # "W_gene_P47Phox_gene_PU1_gene_CEBPa_gene_cRAF"	;	# 136
    100.0 ;  # "W_gene_P47Phox_gene_PPARg"	;	# 137
    100.0 ;  # "W_gene_PPARg_RNAP"	;	# 138
    100.0 ;  # "W_gene_PPARg_gene_Trigger"	;	# 139
    100.0 ;  # "W_gene_PPARg_gene_CEBPa"	;	# 140
    100.0 ;  # "W_gene_PPARg_gene_EGR1"	;	# 141
    100.0 ;  # "W_gene_PPARg_gene_PU1"	;	# 142
    100.0 ;  # "W_gene_PPARg_gene_AP1"	;	# 143
    100.0 ;  # "W_gene_PU1_RNAP"	;	# 144
    100.0 ;  # "W_gene_PU1_gene_Trigger"	;	# 145
    100.0 ;  # "W_gene_PU1_gene_CEBPa"	;	# 146
    100.0 ;  # "W_gene_PU1_gene_PU1"	;	# 147
    100.0 ;  # "W_gene_PU1_gene_AP1"	;	# 148
    100.0 ;  # "W_gene_PU1_gene_OCT1"	;	# 149
    100.0 ;  # "W_gene_PU1_gene_AhR"	;	# 150
    100.0 ;  # "W_gene_PU1_gene_GFI1"	;	# 151
    100.0 ;  # "W_gene_Trigger_RNAP"	;	# 152
    100.0 ;  # "W_gene_cRAF_RNAP"	;	# 153
  ];

  # return the corrected value -
  return parameter_bounds_function(perturbed_parameter_array,lower_bound_array,upper_bound_array);
end

function acceptance_probability_function(rank_array,temperature)
  return (exp(-rank_array[end]/temperature))
end

function cooling_function(temperature)

  # define my new temperature -
  alpha = 0.90
  return alpha*temperature
end


function parameter_bounds_function(parameter_array,lower_bound_array,upper_bound_array)

  # iterate through and fix the parameters -
  new_parameter_array = copy(parameter_array)
  for (index,value) in enumerate(parameter_array)

    lower_bound = lower_bound_array[index]
    upper_bound = upper_bound_array[index]

    if (value<lower_bound)
      new_parameter_array[index] = lower_bound
    elseif (value>upper_bound)
      new_parameter_array[index] = upper_bound
    end
  end

  return new_parameter_array
end
