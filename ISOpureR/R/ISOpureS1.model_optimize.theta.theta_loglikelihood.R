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

### FUNCTION: ISOpureS1.model_optimize.theta.theta_loglikelihood.R #######################################################
#
# Input variables:
#   ww: the theta weights corresponding to patient dd, a 1xK matrix
#   tumordata: a GxD matrix representing gene expression profiles of tumor samples
#   dd: the patient number
#   model: list containing all the parameters to be optimized
# 
# Output variables:
#   loglikelihood: loglikelihood function relevant to optimizing theta

ISOpureS1.model_optimize.theta.theta_loglikelihood <- function(ww, tumordata, dd, model) {
	
	# K = number of normal profiles + 1 
	K <- length(ww); 
	# G = number of genes
	G <- ncol(model$log_all_rates);
   
	# theta_weights for patient dd
	ww <- matrix(ww, nrow=1, ncol=K);
	expww <- exp(ww);
	# theta for patient dd (constrained so that sum of entries is 1) 
	theta <- expww / sum(expww);

	# loglikelihood of theta
	# log p(theta_d | vv) = log Dirichlet( theta_d | vv ), ignoring the vv term in front
	loglikelihood <- (model$vv-1) %*% t(ISOpure.util.matlab_log(theta)); 
 
	# loglikelihood of observed tumour profile t_dd
	# log p(t_d | B, theta_d c_d) = log Multinomial(t_n| alpha*c_d + theta_d*B)
	log_P_t_given_theta <- ISOpure.util.logsum(t(ISOpure.util.repmat(ISOpure.util.matlab_log(theta), G, 1)) + model$log_all_rates,1);
	loglikelihood <- loglikelihood + (log_P_t_given_theta %*% tumordata[,dd]);
	
	# take the negative of the loglikelihood since using a minimizer
	loglikelihood <- -loglikelihood;
	return(as.numeric(loglikelihood));
}
