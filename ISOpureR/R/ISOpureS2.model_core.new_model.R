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

### FUNCTION: ISOpureS2.model_core.new_model.R #############################################################
# 
# Input variables: 
#   tumordata: a GxD matrix representing gene expression profiles of tumor samples
#   kappa: a 1xD matrix which represents strength parameter kappa over cc, given the reference 
#		profile mm
#   INITIAL_VV: a vector with K components, the prior over mixing proportions, theta, with 
#        last entry weighed more heavily 
#   PPtranspose: the prior on the tumor-specific cancer profiles is just the reference cancer
#        profile (1xG matrix) learned in ISOpureS1, standardized so that all entries sum to 1 
#   BBtranspose: a (K-1)xG matrix of the standardized normal profiles, so that they sum to 1
#
# Output variables:
#   model: a newly generated model list to hold all the parameters
#
# REVIST: This does not have backwards compatibility with ISOLATE, as PPtranspose assumed to have 
# just one row. To extend the model to several cancer profiles, may need to modify the way PP
# is treated as a matrix.

ISOpureS2.model_core.new_model <- function(tumordata, kappa, INITIAL_VV, PPtranspose, BBtranspose){

	model <- list();

	# vv vector should be a 1xK matrix
	model$vv <- matrix(INITIAL_VV, nrow=1, ncol=length(INITIAL_VV)); 
	
	model$log_BBtranspose <- as.matrix(log(BBtranspose)); 
	model$PPtranspose <- matrix(PPtranspose, nrow=1, ncol=length(PPtranspose));

	# D = # tumor samples
	# K = # normal profiles + 1 cancer profile
	D <- dim(tumordata)[2];

	model$kappa <- matrix(rep(kappa, D), nrow=1, ncol=D);
	
	# note, in step 2, PPtranspose has just one row
	P <- 1;  
	model$omega <- matrix(1/P, D, P);
	
	# initialize tumor-specific cancer profiles to the reference cancer profile
	model$log_cc <- model$omega %*% PPtranspose;

	return(model)
}
