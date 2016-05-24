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

### FUNCTION: ISOpureS2.model_optimize.theta.theta_loglikelihood.R #######################################################
#
# Input variables:
#   ww: the theta weights correspoding to patient dd, a 1xK matrix
#   tumourdata: a GxD matrix representing gene expression profiles of tumour samples
#   dd: the patient number
#   model: list containing all the parameters to be optimized
# 
# Output variables:
#   loglikelihood: loglikelihood function relevant to optimizing theta

ISOpureS2.model_optimize.theta.theta_loglikelihood <- function(ww, tumordata, dd, model) {
	
	# K = number of normal profiles + 1 
	K <- ncol(model$theta);  
	# G = number of genes
	G <- ncol(model$log_BBtranspose); 

	# theta for patient dd, omitting last column
	# changed the definition of 'remaining' to match Gerald's for numerical stability
	# remaining <- 1-model$theta[dd, K]; 
	remaining <- sum(model$theta[dd,1:(K-1)]);

	ww <- matrix(ww, nrow=1, ncol=length(ww));
	expww <- exp(ww); 
	theta <- remaining * expww / sum(expww);
	
	# all the thetas (including last column which contains tumour purities, i.e.
	# the fraction of cancer cells in the tumour)
	alltheta <- cbind(theta, model$theta[dd,K]);

	# loglikelihood of theta
	# log p(theta_d | vv) = log Dirichlet( theta_d | vv ), ignoring the vv term in front
	loglikelihood <- (model$vv-1) %*% t(ISOpure.util.matlab_log(alltheta)); 
	
	log_all_rates <- rbind(model$log_BBtranspose, model$log_cc[dd,]);

	# loglikelihood of observed tumour profile t_dd
	# log p(t_d | B, theta_d, c_d) = log Multinomial(t_n| alpha*c_d + theta_d*B)

	log_P_t_given_theta <- ISOpure.util.logsum(t(ISOpure.util.repmat(ISOpure.util.matlab_log(alltheta),G, 1)) + log_all_rates,1);

	loglikelihood <- loglikelihood + (log_P_t_given_theta %*% tumordata[,dd]);

	#if (is.finite(loglikelihood)==FALSE) {
	#	stop('something non-finite returned from ISOpure.model_optimize.cg_code.rminimize for theta loglikelihood');
	#}

	# take the negative of the loglikelihood since using a minimizer
	loglikelihood <- -loglikelihood;
	return(as.numeric(loglikelihood));
}
