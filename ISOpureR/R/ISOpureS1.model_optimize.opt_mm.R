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

### FUNCTION: ISOpureS1.model_optimize.opt_mm.R ############################################################################
# 
# Input variables: 
#   tumordata: a GxD matrix representing gene expression profiles of tumour samples
#   model: list containing all the parameters to be optimized
#   NUM_ITERATIONS_RMINIMIZE: minimum number of iteration that the minimization algorithm runs
#   iter: the iteration number
#   NUM_GRID_SEARCH_ITERATIONS: number of times to try restarting with different initial values
#
# Output variables: 
#   model: the model with mm_weights updated
#
# Description: The goal of this function is to optimize the reference cancer profile mm.
#   Because mm is constrained (must be parameters of multinomial/discrete
#   distribution), we don't directly optimize the likelihood function w.r.t.
#   mm, but we perform change of variables to do unconstrained
#   optimization.  We therefore store these unconstrained variables in the
#   field "mm_weights", and update these variables.

ISOpureS1.model_optimize.opt_mm <- function(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS) {
	
	# if mm_weights is not a field (i.e. for the first iteration), initialize
	# mm_weights to the last row of log_all_rates, which in the ISOpure case is just a linear
	# combination of all the normal profiles with equal weights
	if (!any(names(model) == "mm_weights")) {
		model$mm_weights <- model$log_all_rates[nrow(model$log_all_rates),,drop=FALSE];
	}

	init_xx <- t(model$mm_weights);

	retval <- ISOpure.model_optimize.cg_code.rminimize(init_xx, ISOpureS1.model_optimize.mm.mm_loglikelihood, ISOpureS1.model_optimize.mm.mm_deriv_loglikelihood, NUM_ITERATIONS_RMINIMIZE, tumordata=tumordata, model=model);
	
	# mm_weights is the first entry of returned values from ISOpure.model_optimize.cg_code.rminimize, make into a 1xG matrix
	model$mm_weights <- matrix(retval[[1]], 1, length(retval[[1]])); 
	# subtracting ISOpure.util.logsum(mm_weights) from mm_weights corresponds to dividing by total to make sure sum is 1 (once take exponent)
	model$log_all_rates[nrow(model$log_all_rates),] <- model$mm_weights-as.numeric(ISOpure.util.logsum(t(model$mm_weights),1));

	return(model);
}
