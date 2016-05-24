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

### FUNCTION: ISOpureS1.model_optimize.omega.omega_loglikelihood.R #######################################################
#
# Input variables:
#   ww: (K-1)x1 matrix, log(omega), where the entries in omega are constrained to add to 1
#        where K-1 is the number of normal samples
#   tumordata: a GxD matrix representing gene expression profiles of tumor samples
#   model: list containing all the parameters to be optimized
# 
# Output variables:
#   loglikelihood: loglikelihood function relevant to omega

ISOpureS1.model_optimize.omega.omega_loglikelihood <- function(ww, tumordata, model) {
	
	expww <- exp(ww);
	omega <- t(expww) / sum(expww);
	
	kappa <- model$kappa;
	kappaomegaPP <- as.numeric(model$kappa) %*% omega %*% model$PPtranspose;
	
	omegaPP <- omega %*% model$PPtranspose;
	kappaPP <- as.numeric(kappa)  *  model$PPtranspose;

	# log p(m|k,B,omega) = log Dirichlet(m | k*B*omega )
	loglikelihood <- lgamma(sum(kappaomegaPP)) - sum(lgamma(kappaomegaPP)) + ((kappaomegaPP-1) %*% (t(model$log_all_rates[nrow(model$log_all_rates), ,drop=FALSE])));
	
	# take the negative of the loglikelihood since using a minimizer
	loglikelihood <- -loglikelihood;
	return(as.numeric(loglikelihood));
}
