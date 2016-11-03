include("Include.jl")

function updateCholesky(A,zV,CCOV,AVG_P_SUCC,PTHRESH)

  if (AVG_P_SUCC<PTHRESH)
		CA = sqrt(1-CCOV);
		F = 1 + ((1-CA^2)*(norm(zV))^2)/(CA^2);
		A = CA*A+CA/(norm(zV)^2)*(sqrt(F)-1)*A*zV*transpose(zV);
	else
		A = A;
	end

  return A;
end

function updateStepSize(LAMBDA_SUCC,AVG_P_SUCC,P_SUCC_TARGET,CP,D,SIGMA)

	# Compute AVG_P_SUCC -
	AVG_P_SUCC = (1-CP)*AVG_P_SUCC+CP*LAMBDA_SUCC;

	# Compute the step-size -
	F = (1/D)*(AVG_P_SUCC-(P_SUCC_TARGET/(1-P_SUCC_TARGET))*(1-AVG_P_SUCC));
	SIGMA = SIGMA*exp(F);

	if (SIGMA<= 1e-5)
		SIGMA = 1e-5;
	end

  return (SIGMA,AVG_P_SUCC);
end

function sample_function(parameter_guess,kVI,lower_bound_array,upper_bound_array,NPARAMETERS,A,SIGMA)

  # Generate a new vector -
  zV = randn(NPARAMETERS);

  # Compute the new perturbed parameter vector -
  perturbed_parameter_array = parameter_guess.*(1 + SIGMA*A*zV);
  corrected_parameter_array = parameter_bounds_function(perturbed_parameter_array,lower_bound_array,upper_bound_array);

  # return -
  return (corrected_parameter_array,zV)
end

function bounds_function(number_of_parameters)

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
    1000.0 ;  # "K_gene_CEBPa_gene_GFI1"	;	# 22
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
    1000.0 ;  # "K_gene_PU1_gene_PU1"	;	# 82
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

  return (lower_bound_array,upper_bound_array)
end

function estimate_parameters(pObjectiveFunction,initial_parameter_array,data_dictionary,number_of_iterations,error_target)

  # Get the initial parameter set -
  kVI = initial_parameter_array;
  N = length(kVI);

  # Initialize the bounds -
  (lower_bound_array,upper_bound_array) = bounds_function(N);

  # Initialize the algorithm parameters -
  D = 1+N/2;
  P_SUCC_TARGET = 2/11;
  CP = 1/12;
  CC = 2/(N+2);
  CCOV = 2/(N^2+6);
  PTHRESH = 0.99;

  # Set the initial search directions - apply -
  A = eye(N,N);
  AVG_P_SUCC = P_SUCC_TARGET;

  # Setup the run -
  should_loop_continue = true
  parameter_array = kVI;
  SIGMA = 0.10;

  # Initialize -
  error_array = zeros(2)
  error_array[1] = Inf
  counter = 1;
  while (should_loop_continue == true)

    # Generate new parameter set -
  	(perturbed_parameter_array,zV) = sample_function(parameter_array,initial_parameter_array,lower_bound_array,upper_bound_array,N,A,SIGMA);

    # Evaluate the parameter array -
    objective_function_array = pObjectiveFunction(perturbed_parameter_array);
    error_array[2] = sum(objective_function_array);

    # Was the step good or bad?
    delta_error = error_array[2] - error_array[1];
  	if (delta_error<0)
  		LAMBDA_SUCC = 1;
  	else
  		LAMBDA_SUCC = 0;
  	end

  	# Update the step size -
  	(SIGMA,AVG_P_SUCC) = updateStepSize(LAMBDA_SUCC,AVG_P_SUCC,P_SUCC_TARGET,CP,D,SIGMA);

    # Ok, recompute the direction array -
    if (LAMBDA_SUCC == 1)

      # Message -
      msg = "Step improved the objective function. New error = "*string(error_array[2])*" Old error = "*string(error_array[1])*" counter = "*string(counter)
      println(msg)

  		# Keep the peturbed parameters -
  		parameter_array = perturbed_parameter_array;

  		# Update the covariance matrix -
  		A = updateCholesky(A,zV,CCOV,AVG_P_SUCC,PTHRESH);

  		# Update the error -
  		error_array[1] = error_array[2];

  		# If we have a successful step and we are stuck at the min, reset this
  		if (SIGMA<=1e-5)
  			SIGMA = 0.05;
  		end;

  		# Update the counter =
  		counter = counter + 1;

      # Write parameters to disk ...
      writedlm("./parameter_best_v1.dat",parameter_array);

  	else
  		# If we have a successful step and we are stuck at the min, reset this
  		if (SIGMA<=1e-5)
  			SIGMA = 0.05;
  		end;

      msg = "No improvement. Current error = "*string(error_array[1])*" Step error = "*string(error_array[2])*" counter = "*string(counter)
      println(msg)
  	end

  	# Check to see if we need to set FLAG -
  	if (counter>=number_of_iterations || error_array[1] <= error_target)
  		should_loop_continue = false;
  	end
  end

  return parameter_array
end


function main()

  # Initialize -
  number_of_parameters = 153

  # Setup initial parameter guess -
  initial_parameter_array = zeros(number_of_parameters);

  # Load the data dictionary -
  data_dictionary = DataDictionary(0.0,0.0,0.0);

  # Load the parameter array -
  previous_parameter_array = readdlm("./parameter_best_v1.dat.0")

  # Update data dictionary to match new parameters before calculating obj
  parameter_mapping_array = data_dictionary["parameter_name_mapping_array"]
  for index = 1:length(parameter_mapping_array)
    if index <= data_dictionary["number_of_binding"]
      #value = data_dictionary["binding_parameter_dictionary"][parameter_mapping_array[index]]
      value = previous_parameter_array[index]
      initial_parameter_array[index] = value;
    else
      #value = data_dictionary["control_parameter_dictionary"][parameter_mapping_array[index]]
      value = previous_parameter_array[index]
      initial_parameter_array[index] = value;
    end
  end

  # Setup stoping criteria -
  NUMBER_OF_ITERATIONS = 100;
  ERR_TARGET = 10;

  # Setup function pointers -
  pObjectiveFunction = objective_function;

  # estimate parameters -
  parameter_array = estimate_parameters(pObjectiveFunction,initial_parameter_array,data_dictionary,NUMBER_OF_ITERATIONS,ERR_TARGET);

  # Write parameters to disk ...
  writedlm("./parameter_best_v1.dat",parameter_array);
end

# Run me ...
main()
