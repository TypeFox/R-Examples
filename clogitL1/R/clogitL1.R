clogitL1 = function(x, y, strata, numLambda = 100, minLambdaRatio = 0.000001, switch = 0, alpha = 1){
	
	# twiddle parameters for stability
	minLambdaRatio = max(c(minLambdaRatio, 0.000001))
	minLambdaRatio = min(c(minLambdaRatio, 1))

	alpha = max(c(alpha, 0.000001))
	alpha = min(c(alpha, 1))
	
	# make sure inputs are appropriate
	if (length(y) != length(strata)) stop("Please specify one stratum per observation")
	if (length(y) != nrow(x)) stop("Please ensure that each observation has predictors and response")
	if (any(y != 1 & y != 0)) stop("Response vector should contain 1 for cases and 0 for controls")

	# first put y and x into form conducive to call to C function
	# group strata together, putting cases first
	selectOrder = order(strata, 1-y) # select indices such that strata are together, with cases (y=1) coming first
	yC = y[selectOrder]
	xC = as.matrix(x[selectOrder,])
	nVec = table(strata)
	mVec = tapply(y, strata, sum)
	
	# fit the model - clogitL1_c returns a matrix
	# 	first ncol(x) columns: beta estimates
	#	column ncol(x)+1: lambda values, then non zero beta and strong set beta
	esti.beta = clogitl1_c (nVec, mVec, xC, yC, switch, numLambda, minLambdaRatio, alpha=alpha)

	out.beta = esti.beta[,1:ncol(x)] # maybe we stopped before we got to numLambda iterations
	empty.beta = apply (esti.beta, 1, function(x){all(x==0)})
	out.beta = out.beta[!empty.beta,]

	out = list(beta=out.beta, lambda=esti.beta[!empty.beta,ncol(x)+1], nz_beta=esti.beta[!empty.beta,ncol(x)+2], ss_beta=esti.beta[!empty.beta,ncol(x)+3], dev_perc=esti.beta[!empty.beta, ncol(x)+5], y_c=yC, X_c=xC, nVec=nVec, mVec=mVec, alpha=alpha)
	class(out) = "clogitL1"
	out
}
