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

### FUNCTION: ISOpure.util.logsum.R ############################################################################ 
#
# Input variables: 
#   xx: a matrix (not a vector, must be a matrix!)
#   dimen: the dimension along which the long sum is taken 
#         (1 for row, 2 for column)
#
# Output variables:
#    computes ls = log(sum(exp(x),dimen))
#        returns the log of sum of exps, summing over dimension dimen
#        but in a way that tries to avoid underflow/overflow
#
# Notes:  
#   basic idea: shift before exp and reshift back
#   log(sum(exp(x))) = alpha + log(sum(exp(x-alpha)));
#
# REVISIT: comment by Catalina Anghel, 2014-01-17: I asked Gerald Quon why there was the term 
# 2*log(xdims[dimen]) in the computation of alpha, but he didn't know, and said that it could 
# be dropped.

ISOpure.util.logsum <- function(xx, dimen){

	# ensure that the input is a matrix, not a vector
	if (is.null(dim(xx)[dimen])){
		stop('the given dimension of xx is null -- maybe xx is given as a vector instead of a matrix?')
	} 

	if (dim(xx)[dimen]<=1){
		return(xx)
	}

	xdims <- t(as.matrix(dim(xx)));

	## This is in the Matlab code, but it's not correct in Matlab (since function uses the variable 
	## dimen (dim in Matlab) before doing this check.) So eliminate this, and assume that you are 
	## *always* given dimen.
	# if (nargs()<2){
	# 	nonsingletons <- which(xdims>1);
	# 	dimen <- nonsingletons[1];
	# } 
 
	# note command application of max along row/columns is *opposite* of the command in Matlab, ie.
	# max(A,[],2) in Matlab corresponds to apply(A,1,max) in R
	# but since using Rcpp will use Matlab format, ie.
	# if dimen=1, find the  max of each column, if dimen=2 find max of each row
	
	# M is a matrix of the maximum entries of xx or abs(xx) along each row or column
	# Ex: if dimen = 1, then M is the max of each column 
	#     and M = [xx_Max_c1, xx_Max_col2, ..., xx_Max_col_G]
	# need M to be a row matrix
	if (is.complex(xx)) {
		yy <- abs(xx);
		M <- rcppeigen_max_over_columns_or_rows(yy, dimen); 
	}
	else {
		M <- rcppeigen_max_over_columns_or_rows(xx, dimen); 
	}

	# alpha is a modification of M
	# Ex: if dimen = 1, then  alpha = [a_1,  a_2, ... , a_G]
	alpha <- M-log(.Machine$double.xmax)/2+2*log(xdims[dimen]); 
	repdims <- matrix(1,nrow(xdims),ncol(xdims));
	repdims[dimen] <- xdims[dimen];

	# Ex: if dimen = 1, then ISOpure.util.repmat(alpha,repdims[1], repdims[2]) is the matrix
	#     [a_1, a_2, ... , a_G ]
	#     [a_1, a_2, ... , a_G ]
	#     [          ...       ]
	#     [a_1, a_2, ... , a_G ]
	# that's the same size as xx.
	#
	# so exp(xx-ISOpure.util.repmat(alpha,repdims[1], repdims[2])) calculates
	#     [exp(xx_1,1 - a_1), exp(xx_1,2 - a_2), ... , exp(xx_1,G - a_G)]
	#     [exp(xx_2,1 - a_1), exp(xx_2,2 - a_2), ... , exp(xx_2,G - a_G)]
	#     [          ...                                                ]
	#     [exp(xx_K,1 - a_1), exp(xx_K,2 - a_2), ... , exp(xx_K,G - a_G)]
	# 
	# and colSums(exp(xx-ISOpure.util.repmat(alpha,repdims[1], repdims[2])))) is
	#     [ exp(-a_1)( exp(x_1,1) + exp(x_2,1) +  ... + exp(xx_K,1) ),  ... , exp(-a_G) ( exp(x_1,G) + exp(x_2,G) +  ... + exp(xx_K,G) )  ]
	#
	# and ISOpure.util.matlab_log of this is
	#     [ -a_1*log( exp(x_1,1) + exp(x_2,1) +  ... + exp(xx_K,1) ),  ... , -a_G*log( exp(x_1,G) + exp(x_2,G) +  ... + exp(xx_K,G) )  ]
	#
	# so adding alpha to the above gives (finally!) 
	#     [log( exp(x_1,1) + exp(x_2,1) +  ... + exp(xx_K,1),  ... , log( exp(x_1,G) + exp(x_2,G) +  ... + exp(xx_K,G) )  ]
	if (dimen==1){
		logsum <- alpha + ISOpure.util.matlab_log(colSums(exp(xx-ISOpure.util.repmat(alpha,repdims[1], repdims[2]))));
	} else if (dimen==2) {
		logsum <- alpha + ISOpure.util.matlab_log(rowSums(exp(xx-ISOpure.util.repmat(alpha,repdims[1], repdims[2]))));
	}

	return(logsum);
}
