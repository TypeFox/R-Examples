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

### FUNCTION: ISOpureS2.model_optimize.kappa.kappa_loglikelihood.R #######################################################
#
# Input variables:
#   log_kappa: the 1xD matrix log(kappa - model$MIN_KAPPA)
#   model: list containing all the parameters to be optimized
# 
# Output variables:
#   loglikelihood: loglikelihood function relevant to to optimizing kappa
#   
# Notes: CHANGE from step 1: Kappa is a vector

ISOpureS2.model_optimize.kappa.kappa_loglikelihood <- function(log_kappa, model) {

	# recover kappa from log_kappa
	kappa <- (exp(t(log_kappa)) + model$MIN_KAPPA);
	# expww <- exp(log_kappa); 
	
	omegaPP <- model$omega %*% model$PPtranspose; # size (D x G)
	kappaomegaPP <- omegaPP * ISOpure.util.repmat(t(kappa), 1, ncol(model$PPtranspose));  

	# D = number of tumour samples/patients
	D <- nrow(model$log_cc);
	# # G = number of genes - ? not needed 
	# G <- ncol(model$log_cc);

	loglikelihood <- 0;

	# loglikelihood calculation
	for (dd in 1:D) {
		loglikelihood <- loglikelihood + lgamma(sum(kappaomegaPP[dd,])) - sum(lgamma(kappaomegaPP[dd,])) + ((kappaomegaPP[dd,]-1) %*% (model$log_cc[dd,]));
	}

	# take negative because we are using a minimizer
	loglikelihood <- -loglikelihood;
	return(loglikelihood);
}
