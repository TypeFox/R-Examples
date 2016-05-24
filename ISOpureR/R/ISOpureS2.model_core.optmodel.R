# The ISOpureR package is copyright (c) 2014 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

### FUNCTION: ISOpureS2.model_core.optmodel.R #############################################################
# 
# Input variables: 
#   tumordata: a GxD matrix representing gene expression profiles of tumor samples
#   model: list containing all the parameters to be optimized
#   NUM_ITERATIONS (optional): minimum number of iterations of optimization algorithm, default is 35
#
# Output variables: 
#   model: updated model list containing all the parameters
#   This function optimizes the ISOpure parameters cyclically until convergence.

ISOpureS2.model_core.optmodel <- function(tumordata, model, NUM_ITERATIONS=35) {

	total_loglikelihood_old <- -Inf;
	change_ll_frac <- Inf;

	# the number of iterations of the minimization algorithm    
	NUM_ITERATIONS_RMINIMIZE <- 20;

	# for some parameters, may try to run the minimization again 
	# with NUM_GRID_SEARCH_ITERATIONS different initial conditions
	NUM_GRID_SEARCH_ITERATIONS <- 0;
	iter <- 1;

	flog.info("Running ISOpure step 2: PPE - Patient Profile Estimation Step")

	# run for at least 35 iterations
	# if the change in loglikelihood is greater than 1e-7 iterate up to 100 times  
	while (iter <= NUM_ITERATIONS || (change_ll_frac > 0.0000001 && iter < 100)) {
	
		flog.info('--- optimizing cc...');
		model <- ISOpureS2.model_optimize.opt_cc(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS); 

		flog.info('--- optimizing theta...');
		model <- ISOpureS2.model_optimize.opt_theta(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS);

		flog.info('--- optimizing vv...');
		model <- ISOpureS2.model_optimize.opt_vv(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS);

		flog.info('--- optimizing kappa...');
		model <- ISOpureS2.model_optimize.opt_kappa(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS);

		# # the code to optimize omega is here (for future extension), but not
		# # needed because there's only one component to the prior over the
		# # tumor-specific cancer profiles (the reference cancer profile)
		# print('--- optimizing omega...');
		# model <- ISOpureS2.model_optimize.opt_omega(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS);
		
		model$total_loglikelihood <- ISOpureS2.model_core.compute_loglikelihood(tumordata, model);
		change_ll <- model$total_loglikelihood-total_loglikelihood_old;
		change_ll_frac <- abs(change_ll/model$total_loglikelihood);
	
		flog.info('Total log likelihood: %s', model$total_loglikelihood);	
		flog.info('iter: %s/%s, loglikelihood: %s,',  iter, NUM_ITERATIONS, model$total_loglikelihood); 
		flog.info('       change: %s', change_ll_frac);
		
		iter <- iter+1;
		total_loglikelihood_old <- model$total_loglikelihood;
	}

	if (max(abs(model$vv-1)) < 1e-8){
		flog.warn('The values of vv parameter all all near 1')
	}

	return(model)
}
