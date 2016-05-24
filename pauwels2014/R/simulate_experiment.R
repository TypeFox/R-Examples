simulate_experiment <-
function(theta, knobj, experiment_fun){
	## Simulates experiment for parameter theta, with global parameters from knobj
	## and parameter modification and initial condition modified by experiment_fun
	
	theta <- knobj$transform_params(theta)
	temp <- experiment_fun(theta, knobj$global_parameters$initial_conditions)
	
	simulate_experiment_no_transform(temp$theta, temp$initial_conditions, knobj)
}
