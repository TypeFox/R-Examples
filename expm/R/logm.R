### ===== File part of R package expm =====
###
### Function to compute the matrix logarithm
###

logm <- function(x, method = c("Higham08", "Eigen"),
##		 order = 8, trySym = TRUE,
                 tol = .Machine$double.eps)
{
    ## work with "Matrix" too: A<-as.matrix(A)
    d <- dim(x)
    if(length(d) != 2 || d[1] != d[2]) stop("'x' must be a quadratic matrix")

    method <- match.arg(method)
    switch(method,
	   "Higham08" = logm.Higham08(x)
	   ,
	   "Eigen" = {
	       ## AUTHOR: Christophe Dutang
	       ## matrix exponential using eigenvalues / spectral decomposition and
	       ## Ward(1977) algorithm if x is numerically non diagonalisable
	       .Call(do_logm_eigen, x, tol)
	   })
}

