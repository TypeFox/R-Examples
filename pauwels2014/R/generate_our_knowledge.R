generate_our_knowledge <-
function(
		transform_params,
		global_parameters = list(
			atol = 1e-6, # deSolve parameter
			rtol = 1e-6, # deSolve parameter
			tspan = seq(0,100,0.5), # deSolve parameter
			max_step = 50, # deSolve parameter
			max_it = 200, # optimization parameter
			tol = 1e-3, # optimization parameter
			beta = 2, # optimization parameter
			c = 0.0001, # optimization parameter
			n_multi_mod_weight = 20, # sampling parameter
			max_log_like = - 700, # sampling parameter
			centrality_ratio = 0.4, # sampling parameter
			sample_burn_in = 5000, # sampling parameter
			sample_to_keep1 = 10000, # sampling parameter
			sample_step = 1, # sampling parameter
			final_sample = 10000, # sampling parameter
			final_sample_design = 100, # sampling parameter
			n_simu_weights = 100, # risk estimation parameter
			initial_conditions = c(g6 = 1, p6 = 1,p7 = 1,p8 = 1, v6_mrna = 0,v7_mrna = 0,v8_mrna = 0), # Dynamical parameter
			n_params = 9, # Number of free parameters
			param_names =  c("p_degradation_rate", "r6_Kd", "r11_Kd", "pro6_strength", "pro7_strength", "pro9_strength", "rbs6_strength", "rbs7_strength", "rbs8_strength"), # Names of the free parameters
			params = c(p_degradation_rate = 1, r6_Kd = 1, r11_Kd=1, pro6_strength = 1,  pro7_strength = 1, pro9_strength = 1, rbs6_strength = 1, rbs7_strength = 1, rbs8_strength = 1), # An instance of free parameters
			true_params = c(mrna6_degradation_rate =1, mrna7_degradation_rate =1, mrna8_degradation_rate =1, p_degradation_rate = 0.1,  r6_Kd = 2.6, r6_h = 4, r11_Kd=2, r11_h = 2, r12_Kd = 0.2, r12_h = 2, pro6_strength = 1,  pro7_strength = 0.8, pro9_strength = 3.77, rbs6_strength = 5, rbs7_strength = 5, rbs8_strength = 5), # Value of the true parameters
			true_params_T = c(p_degradation_rate = 50, r6_Kd = 56.9162224661803, r11_Kd = 55.0171665943997, pro6_strength = 50, pro7_strength = 48.3848331165323, pro9_strength = 59.6056891700965, rbs6_strength = 61.649500072267, rbs7_strength = 61.649500072267, rbs8_strength = 61.649500072267), ## True parameter value such that transform_params(true_params_T) = true_params (see reverse_params)
			dllname = "pauwels2014"
		),
		datas = list()
	){
  ## A list object to represent what we know about the model
  ##
  ## transform_params allows to switch from log to normal scale and 
  ## also to encode the knowledge we have about some parameters.
  ##
  ## global_parameters lists the simulation and optimization parameters
  ## default value are used in our experiments
  ##
  ## datas gathers measurements that have been carried out as well as 
  ## posterior sample and different risk estimates. 
  
  list(transform_params = transform_params, global_parameters = global_parameters, datas = datas)
  
}
