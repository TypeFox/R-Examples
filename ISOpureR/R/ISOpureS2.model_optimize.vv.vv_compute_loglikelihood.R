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

### FUNCTION: ISOpureS2.model_optimize.vv.vv_compute_loglikelihood.R ##########################################################
#
# Input variables:
#   ww: log(vv-1), a Kx1 matrix
#   sum_log_theta: the column sums of log(theta), a 1xK matrix
#   D: the number of patients (a scalar)
# 
# Output variables:
#   loglikelihood: loglikelihood function relevant to optimizing vv
#
# REVISIT: This function is not used in step 2, so has not been tested

ISOpureS2.model_optimize.vv.vv_compute_loglikelihood <- function(ww, sum_log_theta, D){
	
	# make sure vv is the right shape
	K <- length(ww); # number of normal profiles + 1 
	ww <- matrix(ww, nrow=1, ncol=K);
	vv <- exp(ww) + 1;

	# compute loglikelihood
	loglikelihood <- D*(lgamma(sum(vv)) - sum(lgamma(vv))) +  (vv-1)%*%t(sum_log_theta);
	return(loglikelihood);
}