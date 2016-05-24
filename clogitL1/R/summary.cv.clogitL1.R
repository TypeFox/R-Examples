summary.cv.clogitL1 = function(object, ...){
	# find the minimum CV lambda and coefficient profile
	minInd = which(object$lambda == object$minCV_lambda)
	minCVBeta = object$beta[minInd,]
	minCVNZbeta = object$nz_beta[minInd]

	# find the 1-see-rule lambda and coefficient profile
	minInd = which(object$lambda == object$minCV1se_lambda)
	minCV1seBeta = object$beta[minInd,]
	minCV1seNZbeta = object$nz_beta[minInd]

	list(lambda_minCV=exp(object$minCV_lambda), beta_minCV=minCVBeta, nz_beta_minCV=minCVNZbeta, lambda_minCV1se=exp(object$minCV1se_lambda), beta_minCV1se=minCV1seBeta, nz_beta_minCV1se=minCV1seNZbeta)

}