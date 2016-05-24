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

### FUNCTION: ISOpureS1.model_optimize.kappa.kappa_deriv_loglikelihood.R #########################################################
#
# Input variables:
#   log_kappa: the scalar log(kappa - model$MIN_KAPPA)
#   tumordata: a GxD matrix representing gene expression profiles of tumour samples
#   model: list containing all the parameters to be optimized
# 
# Output variables:
#   deriv_loglikelihood: the derivative of the part of the loglikelihood function relevant to kappa
#   with respect to log kappa (a scalar given that for step 1 of ISOpure kappa is a scalar)
# 
# Description:
#   Instead of performing constrained optimization on kappa directly, we optimize the log of kappa (-MIN_KAPPA?)
#   in an unconstrained fashion.  Thus, if y=log(kappa) and L is the loglikelihood function w.r.t. y,    
#   to optimize L w.r.t. y,   dL/dy = dL/dkappa * dkappa/dy  , where dkappa/dy = exp(y)= exp( log(kappa)).
#   The input into the derivative function is log(kappa - model$MIN_KAPPA).  
# 
# REVISIT: 
#   This function is not vectorized, as FOR STEP 1 in Matlab, as kappa is a scalar.
#   May have to alter this for backwards compatibility with ISOLATE?

ISOpureS1.model_optimize.kappa.kappa_deriv_loglikelihood <- function(log_kappa, tumordata, model) {

	kappa <- as.vector(exp(log_kappa) + model$MIN_KAPPA);
	omega <- as.numeric(model$omega);

	kappaomegaPP <- kappa * t(model$omega) %*% model$PPtranspose; 
	omegaPP <- t(omega) %*% model$PPtranspose;

	# dL/dkappa
	deriv_loglikelihood <- sum(omegaPP)%*%digamma(sum(kappaomegaPP)) - (omegaPP%*%t(digamma(kappaomegaPP))) + model$log_all_rates[nrow(model$log_all_rates), ,drop=FALSE]%*%t(omegaPP);
	
	# dL/dy = dL/dkappa * dkappa/dy where
	# dkappa/dy = exp(y) = exp( log(kappa))
	deriv_loglikelihood <- deriv_loglikelihood * exp(log_kappa)
	
	# take negative of derivative because we are using a minimizer
	deriv_loglikelihood <- -deriv_loglikelihood;
	return(deriv_loglikelihood);
}
