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

### FUNCTION: ISOpureS1.model_core.compute_loglikelihood.R #############################################################
# 
# Input variables: 
#   tumordata: a GxD matrix representing gene expression profiles of tumor samples
#   model: list containing all the parameters updated in ISOpure step one iterations
#
# Output variables: 
#   loglikelihood: the scalar value of the complete loglikelihood obtained given the model parameters

ISOpureS1.model_core.compute_loglikelihood <- function(tumordata, model) {

	loglikelihood <- 0;

	# D is the number of tumor samples
	D <- dim(tumordata)[2];
	# G is the number of transcipts/genes
	G <- dim(model$log_BBtranspose)[2];

	kappaomegaPP <- t(model$omega) %*% model$PPtranspose * model$kappa; 
	# REVISIT: kept kappa as a scalar above, so didn't follow the Matlab code exactly
	
	# loglikelihood of reference cancer profile
	loglikelihood <- loglikelihood + lgamma(sum(kappaomegaPP)) - sum(lgamma(kappaomegaPP)) + ((kappaomegaPP-1) %*% t(model$log_all_rates[nrow(model$log_all_rates), ,drop=FALSE]));

	# loglikelihood of thetas
	loglikelihood <- loglikelihood + D*(lgamma(sum(model$vv))) - sum(lgamma(model$vv)) + sum((model$vv-1) %*% t(ISOpure.util.matlab_log(model$theta)));

	for (dd in 1:D) {
		# loglikelihood of observed tumour profiles t_i
		log_ptgt <- ISOpure.util.logsum(t(ISOpure.util.repmat(ISOpure.util.matlab_log(model$theta[dd,,drop=FALSE]), G, 1)) + model$log_all_rates, 1);
		loglikelihood <- loglikelihood + (log_ptgt %*% tumordata[,dd]);
	}

	return(loglikelihood);
	
}
