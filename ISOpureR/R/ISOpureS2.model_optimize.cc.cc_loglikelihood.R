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

### FUNCTION: ISOpureS2.model_optimize.cc.cc_loglikelihood.R #######################################################
#
# Input variables:
#   ww: the cc_weights for patient dd, with G entries (not sure if a vector or matrix, check)
#   tumordata: a GxD matrix representing gene expression profiles of tumor samples
#   dd: the patient number
#   model: list containing all the parameters to be optimized
# 
# Output variables:
#   loglikelihood: the part of the loglikelihood function relevant to optimizing cc for 
#      patient dd, the cancer profile for that patient

ISOpureS2.model_optimize.cc.cc_loglikelihood <- function(ww, tumordata, dd, model){
	
	# G = number of genes
	G <- length(ww);
	# reshape ww to be a 1xG matrix
	ww <- as.matrix(ww, nrow=1, ncol=G);
	# print(paste('Min of ww is : ', min(ww)));
	log_cancer_rates <- t(ww) - as.numeric(ISOpure.util.logsum(ww,1));

	expww <- t(exp(ww));

	# For ISOpureS2, omega is all 1's model$PPtranspose is mm. 
	kappaomegaPP <- model$kappa[dd] * model$omega[dd,] %*% model$PPtranspose; # 1xG matrix 

	# # K = number of normal samples plus 1
	# K <- ncol(model$theta); 
	# # D = number of patients
	# D <- ncol(tumordata); 
	log_all_rates <- rbind(model$log_BBtranspose, log_cancer_rates); # KxG matrix

	# For patient d,
	# p(c_d| k_d, mm) = Dirichlet(c_d | k_d*mm) = (constant w.r.t. c_d) * prod_(k = 1^G) (c_d)_k ^ ( (k_d)*(mm_k) -1)
	# Hence, if we take the logarithm, the first part is
	# ( (k_d)*(mm_k) -1 ) * log (c_d) 
	loglikelihood <- (kappaomegaPP-1) %*% t(log_cancer_rates);

	# For patient d,
	# p(t_d| B, theta_d, c_d) = Multinomial(t_d | alpha_d*c_d + sum(theta_d_k*B_d_k)) 
	#                         = (constant w.r.t. c_d) * prod_(1^G) (alpha_d*mm + sum(theta_d_k*B_d_k))_ith_component ^(t_d,i)
	# If we take the logarithm, the first part is
	# log (alpha_d*mm + sum(theta_d_k*B_d_k)) 
	# theta contains both theta_d's and alpha_d, and log all rates contains both BB and c_d 

	log_P_t_given_theta <- ISOpure.util.logsum( t(ISOpure.util.repmat( ISOpure.util.matlab_log(t(model$theta[dd,])), G, 1)) + log_all_rates, 1);
	loglikelihood <- as.numeric(loglikelihood) + as.numeric((log_P_t_given_theta %*% tumordata[,dd]));

	# take the negative of the loglikelihood since we're using a minimizer
	loglikelihood <- -loglikelihood;

	if (is.numeric(loglikelihood)==FALSE) {
		stop('imaginary number returned from ISOpure.model_optimize.cg_code.rminimize in ISOpureS2.model_optimize.cc.cc_loglikelihood');
	}

	return(as.numeric(loglikelihood))
}
