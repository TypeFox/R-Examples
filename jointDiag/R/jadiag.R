##
jadiag <- function(M, W_est0=NULL, eps=.Machine$double.eps, itermax=200, 
		keepTrace=FALSE) {
    L <- dim(M)[3]
    d <- dim(M)[1]
    decr <- 1
    logdet <- log(5.184e17)
	W_est <- W_est0
	if (is.null(W_est)) W_est <- diag(d)
    result <- 0
    w <- rep(1,L)
    ctot <- matrix(M,nrow=d)
    ctot.prov <- ctot
    crit <- NA
	myiter <- 0
	B_trace <- W_est
	while (decr > eps & myiter<itermax) {
		if (is.na(logdet))
			stop("log det does not exist")
		res <- .C("jadiagw",as.double(ctot.prov),
			as.double(w),as.integer(d),as.integer(L),as.double(W_est),
			as.double(logdet),as.double(decr),as.double(result))
		ctot.prov <- res[[1]]
		W_est <- res[[5]]
		if(keepTrace) B_trace <- cbind(B_trace,matrix(W_est,nrow=d))
		logdet <- res[[6]]
		decr <- res[[7]]
		crit <- c(crit,res[[8]])
		myiter <- myiter + 1
    }
    W_est <- matrix(W_est,nrow=d)
	if (myiter == itermax) {
		warning("Convergence not reached")
	}
	if (keepTrace) {
		return(list(B=W_est,criter=crit,
			B_trace=array(B_trace,dim=c(d,d,myiter+1))))
	}
	else {
		return(list(B=W_est,criter=crit,B_trace=NULL))
	}
}

