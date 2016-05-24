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

### FUNCTION: ISOpureS2.model_optimize.theta.theta_deriv_loglikelihood.R #########################################################
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

ISOpureS2.model_optimize.theta.theta_deriv_loglikelihood <- function(ww, tumordata, dd, model) {

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
	theta <- remaining*expww / sum(expww);
	
	# all the thetas (including last column which contains tumour purities, i.e.
	# the fraction of cancer cells in the tumour)
	alltheta <- cbind(theta, model$theta[dd,K]);
	
	log_all_rates <- rbind(model$log_BBtranspose, model$log_cc[dd,]);
	log_P_t_given_theta <- ISOpure.util.logsum(t(ISOpure.util.repmat(ISOpure.util.matlab_log(alltheta),G, 1)) + log_all_rates,1);

	# This part is the same as for ISOpureS2.model_optimize.theta.theta_deriv_loglikelihood.R from step 1, except that the last entry in theta is ommitted from the calculation.
	# the first term is d/dtheta log p(theta_d | vv) = log Dirichlet( theta_d | vv ) = d/dtheta (vv-1) log(theta), ignoring the vv term
	# the second term is d/dtheta log p(t_d | B, theta_d, t_d) = d/dtheta log Multinomial (t_d | alpha_d*mm + sum(theta_d_k*B_d_k) )
	dLdtheta <- (model$vv[1:(length(model$vv)-1)]-1)/theta + rowSums(exp(log_all_rates[1:(nrow(log_all_rates)-1),] - ISOpure.util.repmat(log_P_t_given_theta, K-1, 1)) * t(ISOpure.util.repmat(tumordata[,dd],1,K-1)));

	# The dtheta/dw term is similar to the domega/dww term in ISOpureS2.model_optimize.omega.omega_deriv_loglikelihood, where this term is explained in more detail.
	# change of variables
	dthetadw <- -(t(expww)%*%expww)/(sum(expww)^2);

	# code for subdiag from "IShouldBuyABoat" on stackoverflow, accessed Jan 2014: 
	# http://stackoverflow.com/questions/7745363/r-equivalent-to-diagx-k-in-matlab
	subdiag <- function(vec, size, offset=0){ 
	  M <- matrix(0, size, size)
	  M[row(M)-offset == col(M)] <- vec
	  return(M)
	}

	# note offset for subdiag is 0, so the terms expww/sum(expww) are on the diagonal
	dthetadw <- dthetadw + subdiag(expww/sum(expww), length(ww));
	dthetadw <- dthetadw*remaining;
	deriv_loglikelihood <- t(dLdtheta%*%dthetadw);	


	# set the first derivative to be zero, to set the scale of the w's
	deriv_loglikelihood[1] = 0;

	# Rasmussen's conjugate gradient method minimizes, so we take the negative of the derivative
	deriv_loglikelihood <- -deriv_loglikelihood;

	if (is.finite(deriv_loglikelihood)==FALSE) {
		stop('something non-finite returned from ISOpure.model_optimize.cg_code.rminimize in the theta derivative');
	}

	return(as.matrix(deriv_loglikelihood));
}
