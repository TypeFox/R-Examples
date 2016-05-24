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

### FUNCTION: ISOpureS2.model_optimize.opt_kappa.R ############################################################################ 
#
# Input variables: 
#   tumordata: a GxD matrix representing gene expression profiles of tumour samples for each patient
#   model: list containing all the parameters to be optimized
#   NUM_ITERATIONS_RMINIMIZE: minimum number of iteration that the minimization algorithm runs
#   iter: the iteration number
#   NUM_GRID_SEARCH_ITERATIONS: number of times to try restarting with different initial values
#
# Output variables:
#    model: the model with the kappa parameter (which is a 1xD vector) updated
#
# Description:  This function optimizes kappa, the strength parameter in the prior over the reference
# cancer profile.  Note that we don't directly optimize kappa because it has constraints (must be 
# greater than the minimum determined in ISOpure.step2.PPE.)

ISOpureS2.model_optimize.opt_kappa  <- function(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS) {
	
	kappa <- (model$kappa); 
	
	init_xx <- t(log(kappa - model$MIN_KAPPA)); # model$kappa_weights
   
	# perform the optimization
	returnval <- ISOpure.model_optimize.cg_code.rminimize(init_xx, ISOpureS2.model_optimize.kappa.kappa_loglikelihood, ISOpureS2.model_optimize.kappa.kappa_deriv_loglikelihood, NUM_ITERATIONS_RMINIMIZE, model);
	xx <- as.numeric(returnval[[1]]);

	# this makes sure that kappa is at least MIN_KAPPA
	model$kappa <- t(as.matrix(exp(xx) + model$MIN_KAPPA));
	
	return(model);
}
