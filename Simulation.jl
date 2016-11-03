# Runs the model to SS, and issue the pertrubation -
function Simulation(time_start,time_stop,time_step_size,data_dictionary)

  # Run the model to steady-state -
  steady_state_array = EstimateSteadyState(0.01,data_dictionary);

  # Ok, reset the TRIGGER and CRAF -
  steady_state_array[53] = 1;
  steady_state_array[54] = 1;

  #@show steady_state_array[37:54]

  # Reset the IC's -
  data_dictionary["initial_condition_array"] = steady_state_array;

  # Solve the model equations -
  (T,X) = SolveBalances(time_start,time_stop,time_step_size,data_dictionary)

  # return my simulation time, and state -
  return (T,X);
end

function EstimateSteadyState(epsilon,data_dictionary)

  initial_condition_vector = data_dictionary["initial_condition_array"];
  ic_array = copy(data_dictionary["initial_condition_array"])
  number_of_states = length(ic_array)

  # Setup loop -
  EPSILON = epsilon;
  TSTART = 0.0;
  Ts = 1.0;
  TSTOP = 1000;
  did_reach_steady_state = false
  while (!did_reach_steady_state)

    # solve the balances -
    (TSIM,X1) = SolveBalances(TSTART,TSTOP,Ts,data_dictionary)

    # Take a few additional steps -
    TNEXT_START = TSTOP+Ts;
    TNEXT_STOP = TNEXT_START+1.0;
    Ts = 0.1;

    # solve the balances again 0
    initial_condition_array = vec(X1[end,:])
    data_dictionary["initial_condition_array"] = initial_condition_array;
    (TSIM,X2) = SolveBalances(TNEXT_START,TNEXT_STOP,Ts,data_dictionary)

    # Find the difference -
    DIFF = norm((X2[end,:] - X1[end,:]));

    # Should we stop -or- go around again?
    if (DIFF<EPSILON)
      did_reach_steady_state = true;
      return (vec(X2[end,:]));
    else

      # No, we did *not* reach steady state ....
      TSTART = TSTOP+Ts
      TSTOP = 1.0 + TSTART;
      Ts = 0.1;

      initial_condition_array = vec(X2[end,:])
      data_dictionary["initial_condition_array"] = initial_condition_array;
    end
  end

  # return
  return XSS;
end
