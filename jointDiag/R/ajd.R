## a wrapper to joint approximate functions
ajd <- function(M,A0=NULL,B0=NULL,eps=.Machine$double.eps, itermax=200,
		keepTrace=FALSE,methods=c("jedi")) {
	nmeth <- length(methods)
	if (nmeth==1) {
		for (i in 1:nmeth) {
			if (methods=="jedi") 
				res <- jedi(M,A0,eps,itermax,keepTrace)
			if (methods=="uwedge") 
				res <- uwedge(M,B0,eps,itermax,keepTrace)
			if (methods=="jadiag") 
				res <- jadiag(M,B0,eps,itermax,keepTrace)
			if (methods=="ffdiag") 
				res <- ffdiag(M,B0,eps,itermax,keepTrace)
			if (methods=="qdiag") 
				res <- qdiag(M,B0,eps,itermax,keepTrace)
		}
		return(res)
	}

	if (nmeth>1) {
		res <- vector(nmeth,mode="list")
		for (i in 1:nmeth) {
			if (methods[i]=="jedi") 
				res[[i]] <- jedi(M,A0,eps,itermax,keepTrace)
			if (methods[i]=="uwedge") 
				res[[i]] <- uwedge(M,B0,eps,itermax,keepTrace)
			if (methods[i]=="jadiag") 
				res[[i]] <- jadiag(M,B0,eps,itermax,keepTrace)
			if (methods[i]=="ffdiag") 
				res[[i]] <- ffdiag(M,B0,eps,itermax,keepTrace)
			if (methods[i]=="qdiag") 
				res[[i]] <- qdiag(M,B0,eps,itermax,keepTrace)
		}
		names(res) <- methods
		return(res)
	}
}



