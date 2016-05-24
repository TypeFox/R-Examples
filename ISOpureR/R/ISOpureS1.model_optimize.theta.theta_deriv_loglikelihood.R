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

### FUNCTION: ISOpureS1.model_optimize.theta.theta_deriv_loglikelihood.R #########################################################
#
# Input variables:
#   ww: the theta weights correspoding to patient dd, a 1xK matrix
#   tumourdata: a GxD matrix representing gene expression profiles of tumour samples
#   dd: the patient number
#   model: list containing all the parameters to be optimized
# 
# Output variables:
#   deriv_loglikelihood: derivative of the loglikelihood function relevant to optimizing theta, not 
#      with respect to theta but with respect to unconstrained variables

ISOpureS1.model_optimize.theta.theta_deriv_loglikelihood <- function(ww, tumordata, dd, model) {

	# K = number of normal profiles + 1 
	K <- nrow(model$log_all_rates); 
	# G = number of genes
	G <- ncol(model$log_all_rates); 
	
	# theta_weights for patient dd
	ww <- matrix(ww, nrow=1, ncol=K);
	expww <- exp(ww);
	# theta for patient dd (constrained so that sum of entries is 1) 
	theta <- expww / sum(expww);

	# derivative is computed not with respect ot theta directly, but with respect 
	# to unconstrained variables via change of variables
	# this part gives log (alpha_d*mm + sum(theta_d_k*B_d_k))
	log_P_t_given_theta <- ISOpure.util.logsum(t(ISOpure.util.repmat(ISOpure.util.matlab_log(theta),G, 1)) + model$log_all_rates,1); 

	# the first term is d/dtheta log p(theta_d | vv) = log Dirichlet( theta_d | vv ) = d/dtheta (vv-1) log(theta), ignoring the vv term
	# the second term is d/dtheta log p(t_d | B, theta_d, t_d) = d/dtheta log Multinomial (t_d | alpha_d*mm + sum(theta_d_k*B_d_k) )
	dLdtheta <- (model$vv-1)/theta + rowSums(exp(model$log_all_rates - ISOpure.util.repmat(log_P_t_given_theta, K, 1)) * t(ISOpure.util.repmat(tumordata[,dd],1,K)));

	# The following lines calculate dL/dtheta * dtheta/dww (where theta is the same as bb)
	# The code is the same as in ISOpureS1.model_optimize.omega.omega_deriv_loglikelihood for dL/domega * domega/dww.
	# where the second term is explained in more detail. 
	bb <- theta;
	deriv_loglikelihood <- exp(  ISOpure.util.matlab_log(bb) + rep(ISOpure.util.logsum(ISOpure.util.matlab_log(dLdtheta) + ISOpure.util.matlab_log(bb),2), length(bb)) );
	deriv_loglikelihood <- -deriv_loglikelihood + (dLdtheta * bb);
	deriv_loglikelihood <- t(deriv_loglikelihood);

	# fix the first element of the derivative to zero, to fix the scale of the unconstrained variables
	deriv_loglikelihood[1] <- 0;

	# take the negative of the derivative because we are using a minimizer
	deriv_loglikelihood <- -deriv_loglikelihood;
	return(as.matrix(deriv_loglikelihood));
}
