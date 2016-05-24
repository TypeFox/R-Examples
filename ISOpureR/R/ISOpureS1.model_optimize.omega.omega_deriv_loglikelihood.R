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

### FUNCTION: ISOpureS1.model_optimize.omega.omega_deriv_loglikelihood.R ##########################
#
# Input variables:
#   ww: (K-1)x1 matrix, log(omega), where the entries in omega are constrained to add to 1
#        where K-1 is the number of normal samples
#   tumordata: a GxD matrix representing gene expression profiles of tumor samples
#   model: list containing all the parameters to be optimized
# 
# Output variables:
#   deriv_loglikelihood: the derivative of the part of the loglikelihood function relevant to omega
#   with respect to (log) omega 
# 
# Description:
#   Instead of performing constrained optimization on omega directly, we optimize the log of omega
#   in an unconstrained fashion. 

ISOpureS1.model_optimize.omega.omega_deriv_loglikelihood <- function(ww, tumordata, model) {
	
	expww <- exp(ww);
	omega <- t(expww) / sum(expww);

	kappa <- as.numeric(model$kappa);
	kappaomegaPP <- as.numeric(model$kappa) %*% omega %*% model$PPtranspose;

	omegaPP <- omega %*% model$PPtranspose;
	kappaPP <- as.numeric(kappa) * model$PPtranspose; 

	# d ( log p(m|k,B,omega) )/ d(omega) = dL/domega
	dLdomega <- digamma(sum(kappaomegaPP))*as.matrix(rowSums(kappaPP), nrow(kappaPP),1) - (kappaPP%*%t(digamma(kappaomegaPP))) + kappaPP%*%model$log_all_rates[nrow(model$log_all_rates),]

	# Note that dLdomega can contain negative entries... in which case Matlab returns log(z) = log(abs(z)) + 1i*atan2(b,a)
	# thus we use ISOpure.util.matlab_log below.

	# The derivative of the L is dL/domega * domega/dww.
	
	# temp = log sum exp ( log dL/domega + log  omega)
	#      = log sum (dL/domega * omega)
	temp <- ISOpure.util.logsum(t(ISOpure.util.matlab_log(dLdomega)) + ISOpure.util.matlab_log(omega),2);
	# now the exp just gives  [ sum (dL/domega * omega) ] * omega, where * is the dot product
	# sum_{i=1 to K-1} (dL/domega)_i * omega_i) (dot product) omega
	deriv_loglikelihood <- Re(exp(rep(temp,length(omega))+ISOpure.util.matlab_log(omega)));
	# therefore, the derivative becomes
	# - sum_{i=1 to K-1} (dL/domega)_i * omega_i) (dot product) omega + dL/domega (dot product) omega
	deriv_loglikelihood <- -deriv_loglikelihood + (t(dLdomega)*omega);
	deriv_loglikelihood <- t(deriv_loglikelihood);
	# Explanation of the derivative: 
	# Both omega and ww are vectors with K-1 entries.
	# We have: omega = exp ww / sum(exp ww), thus
	# Thus:  omega_i = exp ww_i / sum_{i=1 to K-1} (exp ww_i).
	# And 
	#   d omega_i / d ww_j = omega_i - (omega_i * omega_i), when i = j  
	#                      = - (omega_i*omega_j), when i != j
	# giving a Jacobian domega/dww with the entries on the diagonal which are different from the 
	# the other entires. The expression for the derivative results from the multiplication
	#                                           [ (omega_1 - (omega_1)^2)   (-omega_1 * omega_2)    ...   (-omega_1 * omega_(K-1))       ]
	#                                           [ (-omega_2 * omega_1)      (omega_2 - (omega_2)^2) ...   (-omega_2 * omega_(K-1))       ]
	#                                           [ (-omega_3 * omega_1)      (-omega_3 * omega_2)    ...   (-omega_3 * omega_(K-1))       ]       
	# [dL/domega_1 dL/domega_2 ... dL/domega] * [          ...                                                                           ]
	#                                           [                                                                                        ]
	#                                           [ (-omega_(K-1) * omega_1)  (-omega_(K-1) * omega_1)  ... (-omega_(K-1) - omega_(K-1)^2) ]
	#                                            
		
	# we fix the first element of the derivative to zero, to fix the scale of the 
	# unconstrained variables
	deriv_loglikelihood[1] <- 0;
	
	# take the negative of the derivative since we are using a minimizer
	deriv_loglikelihood <- -deriv_loglikelihood;
	return(as.vector(deriv_loglikelihood));
}
