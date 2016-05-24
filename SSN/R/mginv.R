mginv <-
function(X, tol = sqrt(.Machine$double.eps)) { 
dnx <- dimnames(X) 
if(is.null(dnx)) dnx <- vector("list", 2) 
s <- svd(X) 
nz <- s$d > tol * s$d[1] 
structure( 
if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X, 
dimnames = dnx[2:1]) 
}

