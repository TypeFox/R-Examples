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

### FUNCTION: ISOpureS2.model_optimize.opt_vv.R ############################################################################ 
#
# Input variables: 
#   tumordata: a GxD matrix representing gene expression profiles of tumour samples
#   model: list containing all the parameters to be optimized
#   NUM_ITERATIONS_RMINIMIZE: minimum number of iteration that the minimization algorithm runs
#   iter: the iteration number
#   NUM_GRID_SEARCH_ITERATIONS: number of times to try restarting with different initial values
#
# Output variables:
#    model: the model with the vv parameter updated
#
# Description:  This function optimizes vv, the strength parameter in the prior over the reference
# cancer profile.  Note that we don't directly optimize vv because it has constraints (must be >=1 
# to guarantee real-valued likelihoods).

ISOpureS2.model_optimize.opt_vv <- function(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS) {
	
	# K = number of normal profiles+1
	K <- length(model$vv);  
	# reshape vv into a Kx1 matrix
	vv <- matrix(model$vv, nrow=K, ncol=1)

	# note that we don't directly optimize vvs because it has constraints
	# (must be >=1 to guarantee real-valued likelihoods).
	init_ww <- as.matrix(ISOpure.util.matlab_log(vv-1));
	sum_log_theta <- matrix(colSums(ISOpure.util.matlab_log(model$theta)), nrow=1, ncol=K);

	# perform the optimization
	r_min_ret <- ISOpure.model_optimize.cg_code.rminimize(init_ww, ISOpure.model_optimize.vv.vv_loglikelihood, ISOpure.model_optimize.vv.vv_deriv_loglikelihood, NUM_ITERATIONS_RMINIMIZE, sum_log_theta, ncol(tumordata));
	ww <- matrix(r_min_ret[[1]]);

	# convert back into vv
	vv <- exp(ww)+1;
	model$vv <- matrix(vv,nrow=1, ncol=K);

	return(model);
}
