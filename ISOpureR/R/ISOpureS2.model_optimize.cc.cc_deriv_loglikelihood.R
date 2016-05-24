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

### FUNCTION: ISOpureS2.model_optimize.cc.cc_deriv_loglikelihood.R #######################################################
#
# Input variables:
#   ww: the cc_weights, with G entries (not sure if a vector or matrix, check)
#   tumordata: a GxD matrix representing gene expression profiles of tumor samples
#   dd: the patient number
#   model: list containing all the parameters to be optimized
# 
# Output variables:
#   deriv_loglikelihood: the derivative of the part of the likelihood function relevant to optimizing
#      cc.  The derivative is taken not with respect to vv but with respect to unconstrained variables 
#      via a change of variables

ISOpureS2.model_optimize.cc.cc_deriv_loglikelihood <- function(ww, tumordata, dd, model) {

	# G = number of genes
	G <- length(ww); 

	# reshape ww to be a 1xG matrix
	ww <- as.matrix(ww, nrow=1, ncol=G);
	
	log_cancer_rates <- t(ww) - as.numeric(ISOpure.util.logsum(ww,1));
	expww <- t(exp(ww));
 
	# Q: last multiplication is scalar num x num x (1 x G)?
	kappaomegaPP <- model$kappa[dd] * model$omega[dd,] %*% model$PPtranspose;  

	log_all_rates <- rbind(model$log_BBtranspose, log_cancer_rates); # matrix is KxG

	# derivative is computed not with respect to cc directly, but with respect 
	# to unconstrained variables via change of variables

	# For patient d,
	# p(t_d| B, theta_d, c_d) = Multinomial(t_d | alpha_d*c_d + sum(theta_d_k*B_d_k)) 
	#     = sum( t_d_i )! / prod( t_d_i !) * prod_(1^G) (alpha_d*c_d + sum(theta_d_k*B_d_k))_ith_component ^(t_d,i)
	# The first term in t_d is a constant.  Hence, if we take the logarithm, the first part is

	# log (alpha_d*c_d + sum(theta_d_k*B_d_k))
	# theta contains both theta_d's and alpha_d, and log all rates contains both BB and c_d
	# print(paste('Min of theta is : ', min(model$theta[dd,])));
	# print(paste('Min of long expression is : ', min(min(t(ISOpure.util.repmat( log(t(model$theta[dd,])), G, 1)) + log_all_rates))));

	log_P_t_given_theta <- ISOpure.util.logsum( t(ISOpure.util.repmat( ISOpure.util.matlab_log(t(model$theta[dd,])), G, 1)) + log_all_rates, 1);
	
	# add the derivative of loglikelihood: t_d * 1/(alpha_d*mm + sum(theta_d_k*B_d_k)) * alpha
	# exp(log t_d - log P_t_given_theta + log alpha)
	dLdb <- matrix(0,nrow(log_cancer_rates), ncol(log_cancer_rates)) + Re(exp(ISOpure.util.matlab_log(kappaomegaPP-1)-log_cancer_rates));
	dLdb <- dLdb + exp(t(log(tumordata[,dd])) + log(model$theta[dd,ncol(model$theta)]) - log_P_t_given_theta); # dimension (1 x G)

	# make sure bb sums up to 1
	bb <- expww / sum(expww);
	# The following two lines calculate dL/dbb * dbb/dww. (bb == c_d)
	# The code is the same as in ISOpureS2.model_optimize.omega.omega_deriv_loglikelihood for dL/domega * domega/dww.
	# where the second term is explained in more detail. 

	deriv_loglikelihood <- Re(exp(as.vector(ISOpure.util.logsum(ISOpure.util.matlab_log(dLdb) + ISOpure.util.matlab_log(bb), 2)) + ISOpure.util.matlab_log(bb)));
	deriv_loglikelihood <- -deriv_loglikelihood + (dLdb * bb);

	# fix the first element of the derivative to zero, to fix the scale of the unconstrained variables
	deriv_loglikelihood[1] <- 0;

	# take the negative of the derivative because we are using a minimizer
	deriv_loglikelihood <- t(-deriv_loglikelihood);

	# is.numeric will be false if any entry in deriv_loglikelihood is FALSE
	if (is.numeric(deriv_loglikelihood)==FALSE){   
		stop('imaginary number returned from ISOpure.model_optimize.cg_code.rminimize in ISOpureS2.model_optimize.cc.cc_deriv_loglikelihood');
	}

	return(as.matrix(deriv_loglikelihood));
}
