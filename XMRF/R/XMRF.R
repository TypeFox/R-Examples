XMRF <-
function(X, method="LPGM", stability="bootstrap", N=100, beta=0.01, lmin = 0.01, nlams=20, lambda.path=NULL, parallel=TRUE,nCpus=4, sym=TRUE, th=0.01, sth=0.95, R=max(X), R0=0){

	if(!(method %in% c("LPGM", "PGM", "TPGM", "SPGM", "GGM", "ISM"))){
		ghat = NULL
		return(ghat)
	}
	
	if(method == "LPGM"){
		ghat <- LPGM(X, method="LPGM", stability=stability, N=N, beta=beta, lmin=lmin, nlams=nlams, lambda.path=lambda.path, parallel=parallel, 
						nCpus=nCpus, sym=sym, th=th, sth=sth)
	}
	
	if(method == "PGM"){
		##ghat <- TPGM(X, method="PGM", stability=stability, R=R, N=N, beta=beta, lmin=lmin, nlams=nlams, lambda.path=lambda.path, parallel=parallel, 
		##				ncores=nCpus, sth=sth)
		ghat <- PGM(X, method="PGM", stability=stability, R=R, N=N, beta=beta, lmin=lmin, nlams=nlams, lambda.path=lambda.path, parallel=parallel, 
						ncores=nCpus, sth=sth)
	}
	
	if(method == "TPGM"){
		ghat <- TPGM(X, method="TPGM", stability=stability, R=R, N=N, beta=beta, lmin=lmin, nlams=nlams, lambda.path=lambda.path, parallel=parallel, 
						nCpus=nCpus, sym=sym, th=th, sth=sth)
	}
	
	if(method == "SPGM"){
		ghat <- SPGM(X, method="SPGM", stability=stability, R=R, R0=R0, N=N, beta=beta, lmin=lmin, nlams=nlams, lambda.path=lambda.path, parallel=parallel, 
						nCpus=nCpus, sym=sym, th=th, sth=sth)
	}
	
#	if(method == "WPGM"){
#		ghat <- WPGM(X, method="WPGM", stability=stability, R=R, N=N, beta=beta, lmin=lmin, nlams=nlams, lambda.path=lambda.path, parallel=parallel, 
#						ncores=nCpus, sth=sth)
#	}
	
	if(method == "GGM"){
		ghat <- GGM(X, method="GGM", stability=stability, N=N, beta=beta, lmin=lmin, nlams=nlams, lambda.path=lambda.path, parallel=parallel, 
						nCpus=nCpus, sym=sym, th=th, sth=sth)
	}

	if(method == "ISM"){
		ghat <- ISM(X, method="ISM", stability=stability, N=N, beta=beta, lmin=lmin, nlams=nlams, lambda.path=lambda.path, parallel=parallel, 
						nCpus=nCpus, sym=sym, th=th, sth=sth)
	}
	
	if(!is.null(ghat)){
		ghat$call <- match.call()
	}
	return(ghat)
}
