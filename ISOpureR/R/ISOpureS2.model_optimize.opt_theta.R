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

### FUNCTION: ISOpureS2.model_optimize.opt_theta.R ############################################################################ 
#
# Input variables: 
#   tumordata: a GxD matrix representing gene expression profiles of tumour samples
#   model: list containing all the parameters to be optimized
#   NUM_ITERATIONS_RMINIMIZE: minimum number of iteration that the minimization algorithm runs
#   iter: the iteration number
#   NUM_GRID_SEARCH_ITERATIONS: number of times to try restarting with different initial values
#
# Output variables:
#    model: the model with the theta_weights and theta parameter updated (the first K-1 columns)
#           corresponding to the normal sample contributions

ISOpureS2.model_optimize.opt_theta <- function(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS) {
	
	# K = number of normal samples + 1
	K <- ncol(model$theta);
	# D = number of patients/tumour samples
	D <- nrow(model$theta);
	
	# because thetas are constrained (must be parameters of multinomial/discrete
	# distribution), we don't directly optimize the likelihood function w.r.t.
	# theta, but we perform change of variables to do unconstrained
	# optimization.  We therefore store these unconstrained variables in the
	# field "theta_weights", and update these variables
	
	# Furthermore, note that we are fixing the tumor purities (last column of
	# theta), so we are only storing/updating the remaining columns of theta,
	# and optimizing them to sum to 1-alpha_i for tumor i
	if (!any(names(model)=='theta_weights')) {
		model$theta_weights <- log(model$theta[ ,1:(K-1)]);
		}
	
	# update each theta_d separately
	for (dd in 1:D) {
		init_xx <- t(model$theta_weights[dd, ,drop=F]);

		# changed the definition of 'remaining' to match Gerald's change for numerical stability
		# remaining <- 1-model$theta[dd, K]; 
		remaining <- sum(model$theta[dd,1:(K-1)]);

		# perform the optimization
		returnval <- ISOpure.model_optimize.cg_code.rminimize(init_xx, ISOpureS2.model_optimize.theta.theta_loglikelihood, ISOpureS2.model_optimize.theta.theta_deriv_loglikelihood, NUM_ITERATIONS_RMINIMIZE, tumordata=tumordata,dd=dd,model=model)
		xx <- returnval[[1]]; 

		# convert from unconstrained variables to theta
		# constrain theta to sum up to 1
		model$theta[dd, 1:(K-1)] <- remaining * (t(matrix(exp(xx)))/sum(exp(xx)));
		model$theta_weights[dd,] = t(matrix(xx));
	}
	
	return(model);
}
