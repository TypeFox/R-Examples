#-----------------------------------------------------------------------#
# Package: sparc      	   		   									    #
# Method: Semiparametric Generalized Linear Models				 		#
# Authors: Tuo Zhao and Han Liu            								#
# Date: Aug 31st 2013                                                   #
# Version: 0.9.0 	                                                  	#
#-----------------------------------------------------------------------#

sparc = function(X, y, lambda = NULL, lambda.min.ratio=NULL, nlambda = NULL, thol = 1e-4, max.ite = 1e4, alpha = sqrt(1/2)){
	
	n = nrow(X)
	d = ncol(X)
	
	N = n*(n-1)/2
	
	out.SGLM = .C("SGLM", A = as.double(rep(0,N*d)), X = as.double(X), y = as.double(y), NN = as.integer(N), nn = as.integer(n), dd = as.integer(d),lambda = as.double(rep(0,d)),package="sparc")
	
	A = matrix(out.SGLM$A,ncol=d)
	
	if(is.null(lambda)){
		nlambda = 30
		lambda.min.ratio = 0.1
		lambda.max = max(abs(out.SGLM$lambda))
		lambda = exp(seq(0,log(lambda.min.ratio),length=nlambda+1))*lambda.max
		lambda = lambda[-1]
	} else {
		nlambda = length(lambda)
	}
	
	L0 = max(eigen(t(A)%*%A,only.values=T)$values)
	
	out.SLRM = .C("SLRM", A = as.double(A), lambda = as.double(lambda), nnlambda = as.integer(nlambda),LL0 = as.double(L0), nn = as.integer(N), dd = as.integer(d), ww = as.double(rep(0,d*nlambda)), mmax_ite = as.integer(max.ite), tthol = as.double(thol), aalpha = as.double(alpha),package="sparc")
	
	w = matrix(out.SLRM$ww,ncol=nlambda)
	
	return(w)
}