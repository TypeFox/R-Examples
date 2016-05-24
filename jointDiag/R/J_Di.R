##
jedi <- function(M,A0=NULL,eps=.Machine$double.eps, itermax=200,
		keepTrace=FALSE) {
	mat_size <- dim(M)[1]
	nb_mat <- dim(M)[3]
	A <- A0
	if (is.null(A)) {
		A <- diag(mat_size)
	}
	A_trace <- A
	iter <- 0
	crit <- NA
	sh_max <- s_max <- Inf
	while ((iter<itermax) & ((sh_max>eps) | (s_max>eps))) {
		res <- .C("sweepjedi",as.double(M),as.integer(mat_size),
			as.integer(nb_mat),as.double(0),as.double(0),as.double(A))
		A <- array(res[[6]],dim=c(mat_size,mat_size))
		M <- array(res[[1]],dim=c(mat_size,mat_size,nb_mat))
		s_max <- res[[4]]
		sh_max <- res[[5]]
		iter=iter+1
		crit <- c(crit,max(s_max,sh_max))
		if(keepTrace) A_trace <- cbind(A_trace,A)
	}
	if (iter == itermax) {
		warning("Convergence not reached")
	}
	if (keepTrace) {
		return(list(A=A,criter=crit,
			A_trace=array(A_trace,dim=c(mat_size,mat_size,iter+1))))
	}
	else {
		return(list(A=A,criter=crit,A_trace=NULL))
	}
}







