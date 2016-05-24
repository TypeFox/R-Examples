ISM <-
function(X,method="ISM",stability="bootstrap", N=100,beta=0.05, lmin = 0.01, nlams=20, lambda.path=NULL ,parallel=T,nCpus=4,sym=TRUE,th=0,sth=0.8){
	ghat <- glm.generic(X, method=method, link="binomial", stability=stability, N=N, beta=beta, lmin=lmin, nlams=nlams, lambda.path=lambda.path, parallel=parallel, nCpus=nCpus,sym=sym,th=th,sth=sth)
	if(!is.null(ghat)){
		ghat$call <- match.call()
	}
	
	return(ghat)
}
