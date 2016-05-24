simulate_experiment_no_transform <-
function(theta, initial_conditions, knobj){
	## Simulates experiment for parameter theta, with global parameters from knobj.
	## No transformation is applied. This is just a call to DeSolve.
	
	ode(knobj$global_parameters$initial_conditions, method = "ode23", knobj$global_parameters$tspan , func = "derivs", parms = theta, dllname = knobj$global_parameters$dllname, initfunc = "initmod", nout = 1, outnames = "out", rtol = knobj$global_parameters$rtol, atol = knobj$global_parameters$atol, maxsteps = knobj$global_parameters$max_step)
}
