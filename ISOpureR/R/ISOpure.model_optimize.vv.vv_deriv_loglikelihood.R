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

### FUNCTION: ISOpure.model_optimize.vv.vv_deriv_loglikelihood.R ############################################################
#
# Input variables:
#   ww: log(vv-1), a Kx1 matrix
#   sum_log_theta: the column sums of log(theta), a 1xK matrix
#   DD: the number of patients (a scalar)
# 
# Output variables:
#   deriv_loglikelihood: derivative of the loglikelihood function relevant to vv

ISOpure.model_optimize.vv.vv_deriv_loglikelihood <- function(ww, sum_log_theta, DD) {
	
	# K = number of normal profiles + 1 
	K <- length(ww); 
	ww <- matrix(ww, nrow=1, ncol=K);
	sum_log_theta <- matrix(sum_log_theta, nrow=1, ncol=K);
	# a is actually the current vv vector
	a <- exp(ww) + 1;

	# d ( log Dirichlet(theta_n|vv) ) / d(vv) * d(vv)/d(ww), same as for the kappa
	deriv_loglikelihood <- (DD * (digamma(sum(a)) - digamma(a)) + sum_log_theta) * exp(ww);
	# take negative because we are using a minimizer
	deriv_loglikelihood <- t(-deriv_loglikelihood);

	return(as.vector(deriv_loglikelihood));
}
