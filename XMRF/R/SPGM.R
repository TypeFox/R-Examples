SPGM <-
function(X, method="SPGM", stability="bootstrap", R, R0=0, N=100,beta=0.05, lmin = 0.01, nlams=20, lambda.path=NULL ,parallel=TRUE,nCpus=4, sym=TRUE,th=0, sth=0.8){
	if (R < 0){
		cat("ERROR: Truncating threshold R should be positive. \n")
		ghat = NULL
		return(ghat)
	}
	
	# Generate the matrix with values sublinearly truncated between R (upepr bound) and R0 (lower bound)
	Xorig <- X
	X <- round(Bsublin(X, R, R0))
	
	ghat <- LPGM(X, method="SPGM", stability="bootstrap", N=N, beta=beta, 
					lmin=lmin, nlams=nlams, lambda.path=lambda.path, parallel=parallel, nCpus=nCpus,sym=sym,th=th,sth=sth)

	return(ghat)
}
