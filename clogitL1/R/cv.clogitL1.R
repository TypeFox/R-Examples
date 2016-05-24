cv.clogitL1 = function(clObj, numFolds=10){
	# Function that accepts a object of type "clogitL1" and returns an object 
	# of type "cv.clogitL1" containing relevant CV information

	if (class(clObj) != "clogitL1"){
		stop("requires object of type clogitL1")
	}
		
	# assign folds
	K = length(clObj$nVec) # number of folds
	nFolds = min(c(K, numFolds))
	
	perFold = floor(K/nFolds)
	extra = K %% nFolds
	foldVec = NULL
	for (i in 1:nFolds){
		foldVec = c(foldVec, rep(i, perFold))
		if (i <= extra) foldVec = c(foldVec, i) # another one to make up extra
	}
	foldVec = sample (foldVec, K, replace = F) # jumble around the folds

	# now run the CV version of clogitL1
	cv.mat = NULL
	for (i in 1:nFolds){
		keepvec = ifelse (foldVec == i, 0, 1)
		esti.beta.cv = clogitl1CV_c (clObj$nVec, clObj$mVec, clObj$X_c, clObj$y_c, clObj$lambda, keepvec, clObj$alpha)
		cv.mat = cbind(cv.mat, esti.beta.cv[,ncol(clObj$X_c) + 6])
	}

	# find minimum and 1-se rule
	mean.cv.err = -apply (cv.mat, 1, mean)[-1] # cv error at each lambda
	sd.cv.err = apply (cv.mat, 1, function(x){sqrt(var(x))})[-1]/sqrt(10) # estimate of standard error of this mean cv error
	plot.lambdas = log(clObj$lambda)[-1]
	
	min.lambda = plot.lambdas[which.min(mean.cv.err)] # lambda of minimum cv error
	min.lambda.1sd = max(plot.lambdas[mean.cv.err-sd.cv.err <= min(mean.cv.err)])

	out = list(cv_dev=cv.mat, lambda=plot.lambdas, folds=foldVec, mean_cv=mean.cv.err, se_cv=sd.cv.err, minCV_lambda=min.lambda, minCV1se_lambda=min.lambda.1sd, nz_beta=clObj$nz_beta[-1], beta=clObj$beta[-1,])
	class(out) = "cv.clogitL1"
	out
}

