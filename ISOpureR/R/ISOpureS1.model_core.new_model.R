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

### FUNCTION: ISOpureS1.model_core.new_model.R #############################################################
# 
# Input variables: 
#   tumordata: a GxD matrix representing gene expression profiles of tumor samples
#   kappa: scalar strength parameter kappa placed over the reference cancer profile m
#   INITIAL_VV: a vector with K components, the prior over mixing proportions, theta, with last 
# 		entry weighed more heavily 
#   PPtranspose: a (K-1)xG matrix, standardized so that all entries sum to 1, see 
#		ISOpure.step1.CPE.R 
#   BBtranspose: a (K-1)xG matrix of the standardized normal profiles, so that they sum to 1
#
# Output variables: 
#   model: a newly generated model list to hold all the parameters vv, log_BBtranspose, PPtranspose,
#          kappa, theta, omega, log_all_rates

ISOpureS1.model_core.new_model <- function(tumordata, kappa, INITIAL_VV, PPtranspose, BBtranspose){

	model <- list();

	# vv vector should be a 1xK matrix; however using as.matrix() makes the dimension Kx1
	# hence the additional transposition
	model$vv <- t(as.matrix(INITIAL_VV)); 
	
	model$log_BBtranspose <- as.matrix(log(BBtranspose)); 
	model$PPtranspose <- as.matrix(PPtranspose);
	model$kappa <- kappa;

	# D = # tumor samples
	# K = # normal profiles + 1
	D <- dim(tumordata)[2];
	K <- length(model$vv);

	# randomly initialize theta but give higher weight to cancer component initially
	model$theta <- matrix(0, D, K);
	model$theta[,ncol(model$theta)] <- 0.5;
	model$theta[,1:(ncol(model$theta)-1)] <- 0.5/(K-1);
	model$theta <- as.matrix(model$theta / ISOpure.util.repmat(rowSums(model$theta), 1, K)); 

	# P should be K-1, number of normal profiles
	P <- nrow(PPtranspose);
	model$omega <- matrix(1/P, P, 1);
	
	# model.log_all_rates = [B mm]^T , where B is the matrix of normal profiles
	# [b_1 ... b_(K-1)] , and mm is the reference cancer profile
	model$log_all_rates <- as.matrix(rbind(model$log_BBtranspose, log(t(model$omega)%*%model$PPtranspose)));

	return(model)
}
