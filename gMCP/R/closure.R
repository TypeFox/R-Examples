fastgMCP <- function(m, w, p, a, keepWeights=TRUE) {
	if (length(a)>1) {
		warning("Only the first value of 'a' is used!")
	}
	n <- dim(m)[1]
	if (dim(m)[2]!=n || length(w)!=n || length(p)!=n) {
		stop("Wrong dimensions in fastgMCP call!")
	}	
	result <- .C("cgMCP", oldM=as.double(m), oldW=as.double(w), 
			p=as.double(p), a=as.double(a), n=as.integer(n), 
			s=double(n), m = double(n*n), w = double(n))
	return(list(m=matrix(result$m, nrow=n), w=result$w, rejected=(result$s==1)))	
}
