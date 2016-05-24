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

### FUNCTION: ISOpureS1.model_optimize.kappa.kappa_loglikelihood.R #######################################################
#
# Input variables:
#   log_kappa: the scalar log(kappa - model$MIN_KAPPA)
#   tumourdata: a GxD matrix representing gene expression profiles of tumour samples
#   model: list containing all the parameters to be optimized
# 
# Output variables:
#   loglikelihood: loglikelihood function relevant to kappa
# 
# REVISIT: 
#   This function is not vectorized, as FOR STEP 1 in Matlab, as kappa is a scalar.
#   May have to alter this for backwards compatibility with ISOLATE?


ISOpureS1.model_optimize.kappa.kappa_loglikelihood <- function(log_kappa, tumordata, model) {

	# recover kappa from log_kappa
	# note: need as.numeric (or as.vector) here
	# otherwise the matrix multiplication for kappaomegaPP won't work 
	kappa <- as.numeric(exp(log_kappa) + model$MIN_KAPPA); 
	
	omega <- as.numeric(model$omega);

	kappaomegaPP <- as.numeric(kappa * t(model$omega) %*% model$PPtranspose); 
	
	# loglikelihood calculation
	loglikelihood <- lgamma(sum(kappaomegaPP)) - sum(lgamma(kappaomegaPP)) + ((kappaomegaPP-1) %*% (t(model$log_all_rates[nrow(model$log_all_rates), ,drop=FALSE])));

	# take negative because we are using a minimizer
	loglikelihood <- -loglikelihood;
	return(loglikelihood);
}
