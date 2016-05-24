##' @description Generate a random correlation matrix from specified eigenvalues. If eigenvalues are
##' not specified, they are randomly generated from a uniform [0,10] distribution.  
##'
##' @title Generate a random correlation matrix
##' @param n the number of rows/dimensions of the correlation matrix 
##' @param ev Eigenvalues. Defaults to sampling from a uniform distribution betwen 0 and 10. 
##' @return a correlation matrix, of size \code{nxn}
##' @export
##' @references https://stat.ethz.ch/pipermail/r-help/2008-February/153708.html
##' @author Dustin Fife
random.correlation = function(n, ev=runif(n,0,10)){
	Z = matrix(ncol=n, rnorm(n^2))
	decomp = qr(Z)
	Q = qr.Q(decomp) 
	R = qr.R(decomp)
	d = diag(R)
	ph = d / abs(d)
	O = Q %*% diag(ph)
	Z = t(O) %*% diag(ev) %*% O
	return(cov2cor(Z))
}
