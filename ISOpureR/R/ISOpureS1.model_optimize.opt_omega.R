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

### FUNCTION: ISOpureS1.model_optimize.opt_omega.R ############################################################################ 
#
# Input variables: 
#   tumordata: a GxD matrix representing gene expression profiles of tumour samples
#   model: list containing all the parameters to be optimized
#   NUM_ITERATIONS_RMINIMIZE: minimum number of iteration that the minimization algorithm runs
#   iter: the iteration number
#   NUM_GRID_SEARCH_ITERATIONS: number of times to try restarting with different initial values
#
# Output variables:
#    model: the model with the omega_weights and omega parameters updated
#
# Description:  This function optimizes omega, in fact the convex mixing weights that govern prior 
# over the reference cancer profile.
#
# REVISIT: Re-test NUM_GRID_SEARCH_ITERATIONS > 0 part, as had some issues with it in initial testing. 

ISOpureS1.model_optimize.opt_omega  <- function(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS) {
	
	# P is the number of normals in this "site of origin panel"
	P <- length(model$omega);

	# because omega variables are constrained (must be parameters of multinomial/discrete
	# distribution), we don't directly optimize the loglikelihood function w.r.t.
	# omega, but we perform change of variables to do unconstrained
	# optimization.  We therefore store these unconstrained variables in the
	# field "omega_weights", and update these variables
	if (!any(names(model)=='omega_weights')) {
		model$omega_weights <- log(model$omega);
	}   

	init_xx <- model$omega_weights;
	
	# perform the optimization
	returnval <- ISOpure.model_optimize.cg_code.rminimize(init_xx, ISOpureS1.model_optimize.omega.omega_loglikelihood, ISOpureS1.model_optimize.omega.omega_deriv_loglikelihood, NUM_ITERATIONS_RMINIMIZE, tumordata=tumordata, model=model);
	xx <- matrix(returnval[[1]],nrow=P,ncol=1);

	# convert from unconstrained variables to omega
	model$omega <- exp(xx)/sum(exp(xx));
	model$omega_weights <- xx;
	loglikelihood <- ISOpureS1.model_optimize.omega.omega_compute_loglikelihood(model$omega, tumordata, model);
	
	# do some random restarts for omega_weights optimization  
	if (iter <= NUM_GRID_SEARCH_ITERATIONS) {
	
		 for (ii in seq(0,10,1)) {
			# add 0.01 to make sure no number is 0 (or really close to 0)
			init_xx <- matrix(runif(P) + 0.01, nrow=P, ncol=1); 

			returnval <- ISOpure.model_optimize.cg_code.rminimize(log(init_xx), ISOpureS1.model_optimize.omega.omega_loglikelihood, ISOpureS1.model_optimize.omega.omega_deriv_loglikelihood, NUM_ITERATIONS_RMINIMIZE, tumordata=tumordata, model=model);
			xx <- matrix(returnval[[1]],nrow=P,ncol=1);
			newomega <- exp(xx)/sum(exp(xx));
			newll <- ISOpureS1.model_optimize.omega.omega_compute_loglikelihood(newomega, tumordata, model);
		
			if (newll > loglikelihood) {
				model$omega_weights <- xx;
				model$omega <- newomega;
				loglikelihood <- newll;
			}
		}
	}
	
	return(model);
}
