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

### FUNCTION: ISOpureS1.model_optimize.mm.mm_loglikelihood.R #######################################################
#
# Input variables:
#   ww: the mm_weights, with G entries (not sure if a vector or matrix, check)
#   tumordata: a GxD matrix representing gene expression profiles of tumour samples
#   model: list containing all the parameters to be optimized
# 
# Output variables:
#   loglikelihood: the part of the loglikelihood function relevant to optimizing mm, the 
#       reference cancer profile

ISOpureS1.model_optimize.mm.mm_loglikelihood <- function(ww, tumordata, model) {
	
	# # K = number of normal profiles + 1
	# K <- dim(model$log_all_rates)[1]
	# G = number of genes
	G <- dim(model$log_all_rates)[2] 
	# D = number of patients
	D <- ncol(tumordata);
	
	# reshape ww to be a 1xG matrix
	ww <- as.matrix(ww, nrow=1, ncol=G);
	
	log_cancer_rates <- t(ww) - as.numeric(ISOpure.util.logsum(ww,1));
	expww <- exp(t(ww));
	
	kappaomegaPP <- as.numeric(model$kappa) %*% t(model$omega) %*% model$PPtranspose;
	log_all_rates <- rbind(model$log_BBtranspose, log_cancer_rates);
	loglikelihood <- (kappaomegaPP-1) %*% t(log_cancer_rates);

	for (dd in 1:D) {
		# For patient d,
		# p(t_d| B, theta_d, mm) = Multinomial(t_d | alpha_d*mm + sum(theta_d_k*B_d_k)) 
		#     = sum( t_d_i )! / prod( t_d_i !) * prod_(1^G) (alpha_d*mm + sum(theta_d_k*B_d_k))_ith_component ^(t_d,i)
		# The first term in t_d is a constant w.r.t. mm.  Hence, if we take the logarithm, the first part is

		# log (alpha_d*mm + sum(theta_d_k*B_d_k)) 
		# theta contains both theta_d's and alpha_d, and log all rates contains both BB and mm  
		log_P_t_given_theta <- ISOpure.util.logsum( t(ISOpure.util.repmat( ISOpure.util.matlab_log(t(model$theta[dd,])), G, 1)) + log_all_rates, 1);

		# and then we multiply this (dot product) with t_d and add it to the previous loglikelihood
		loglikelihood <- as.numeric(loglikelihood) + as.numeric((log_P_t_given_theta %*% tumordata[,dd]));
	}
	# at the end of the loop, we get a sum over all the patients dd

	# take the negative of the loglikelihood since we're using a minimizer
	loglikelihood <- -loglikelihood;  
	return(as.numeric(loglikelihood));
}
