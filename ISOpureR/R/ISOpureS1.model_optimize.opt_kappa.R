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

### FUNCTION: ISOpureS1.model_optimize.opt_kappa.R ############################################################################ 
#
# Input variables: 
#   tumourdata: a GxD matrix representing gene expression profiles of tumour samples
#   model: list containing all the parameters to be optimized
#   NUM_ITERATIONS_RMINIMIZE: minimum number of iteration that the minimization algorithm runs
#   iter: the iteration number
#   NUM_GRID_SEARCH_ITERATIONS: number of times to try restarting with different initial values
#
# Output variables:
#    model: the model with the kappa parameter updated
#
# Description:  This function optimizes kappa, the strength parameter in the prior over the reference
# cancer profile.  Note that we don't directly optimize kappa because it has constraints (must be 
# greater than the minimum determined in ISOpure.step1.CPE.)
#
# REVISIT: 
# 1. This function is not vectorized, as FOR STEP 1 in Matlab, as kappa is a scalar.
#    May have to alter this for backwards compatibility with ISOLATE?
# 2. Re-test NUM_GRID_SEARCH_ITERATIONS > 0 part, as had some issues with it in initial testing. 

ISOpureS1.model_optimize.opt_kappa  <- function(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS) {
	
	kappa <- as.numeric(model$kappa); 
	
	# note, kappa in our case is a scalar, so I didn't reproduce the matrix operations (transpose, etc) from
	# the Matlab code in this R code
	init_xx <- log(kappa - model$MIN_KAPPA); # model$kappa_weights
	
	# perform the optimization
	returnval <- ISOpure.model_optimize.cg_code.rminimize(init_xx, ISOpureS1.model_optimize.kappa.kappa_loglikelihood, ISOpureS1.model_optimize.kappa.kappa_deriv_loglikelihood, NUM_ITERATIONS_RMINIMIZE, tumordata=tumordata, model=model);
	xx <- as.numeric(returnval[[1]]);

	# this makes sure that kappa is at least MIN_KAPPA
	model$kappa <- exp(xx) + model$MIN_KAPPA;  
	
	# do some random restarts for kappa optimization  
	if (iter <= NUM_GRID_SEARCH_ITERATIONS) {
	    loglikelihood <- ISOpureS1.model_optimize.kappa.kappa_compute_loglikelihood(model$kappa, tumordata, model);

	    MIN_POW_KAPPA <- ceiling(log10(model$MIN_KAPPA));
	
	    for (pow in MIN_POW_KAPPA:15) {
	        # again, this part is vectorized in the Matlab
	        scales <- 10^pow * rep(1, length(model$kappa));
	        init_xx <- log(scales - model$MIN_KAPPA);
	
	        returnval <- ISOpure.model_optimize.cg_code.rminimize(init_xx, ISOpureS1.model_optimize.kappa.kappa_loglikelihood, ISOpureS1.model_optimize.kappa.kappa_deriv_loglikelihood, NUM_ITERATIONS_RMINIMIZE, tumordata=tumordata, model=model);
	        newkappa <- exp(returnval[[1]]) + model$MIN_KAPPA;
	        newll <- ISOpureS1.model_optimize.kappa.kappa_compute_loglikelihood(newkappa, tumordata, model);

	        if (newll > loglikelihood) {
	            model$kappa <- newkappa;
	            loglikelihood <- newll;
	        }
	    }
	}
	
	return(model);
}
