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

### FUNCTION: ISOpureS2.model_optimize.kappa.kappa_deriv_loglikelihood.R #########################################################
#
# Input variables:
#   log_kappa: the 1xD matrix log(kappa - model$MIN_KAPPA)
#   model: list containing all the parameters to be optimized
# 
# Output variables:
#   deriv_loglikelihood: the derivative of the part of the loglikelihood function relevant to kappa
#   with respect to log kappa 

ISOpureS2.model_optimize.kappa.kappa_deriv_loglikelihood <- function(log_kappa, model) {

	kappa <- (exp(t(log_kappa)) + model$MIN_KAPPA);
	expww <- exp(log_kappa); 

	omegaPP <- model$omega %*% model$PPtranspose;
	kappaomegaPP <- omegaPP * ISOpure.util.repmat(t(kappa), 1, ncol(model$PPtranspose));  # (?) the kappa is transposed in the Matlab code

	D <- nrow(model$log_cc);
	G <- ncol(model$log_cc);

	deriv_loglikelihood <- 0; 
	
	# dL/dkappa
	for (dd in 1:D) {
		deriv_loglikelihood[dd] <-  digamma(kappa[dd]) - (omegaPP[dd,,drop=FALSE]%*%t(digamma(kappaomegaPP[dd,, drop=FALSE]))) +( omegaPP[dd,,drop=FALSE]%*%t(model$log_cc[dd,,drop=FALSE]));
	}

	# dL/dy = dL/dkappa * dkappa/dy where
	# dkappa/dy = exp(y) = exp( log(kappa))
	deriv_loglikelihood <- deriv_loglikelihood * expww;
	
	# take negative of derivative because we are using a minimizer
	deriv_loglikelihood <- -deriv_loglikelihood;
	return(as.matrix(deriv_loglikelihood, nrow=D, ncol=1));
}
