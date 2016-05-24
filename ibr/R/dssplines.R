DuchonQ <- function(x,xk,m=2,s=0,symmetric=TRUE) {
  p <- ncol(x)
  nx <- nrow(x)
  if (!symmetric) {
    if (ncol(xk)!=p) stop("number of variables of x and xk must be the same\n")
    nxk <- nrow(xk)
  } else {
    nxk <- nx
    xk <- 0
  }
  k <- 2*m + 2*s - p
  negatif <- if ((1-2*((floor(k/2)+1)%%2))==-1) 1 else 0
  if (k%%2==0) {
    res <- .C("semikerlog",as.double(x),as.double(xk),as.integer(nx),as.integer(nxk),as.double(k/2),as.integer(p),as.integer(negatif),E=double(nx*nxk),as.integer(symmetric))
   } else {
    res <- .C("semikerpow",as.double(x),as.double(xk),as.integer(nx),as.integer(nrow(xk)),as.double(k/2),as.integer(p),as.integer(negatif),E=double(nx*nxk),as.integer(symmetric))
  }
 return(matrix(res$E,nrow=nx,ncol=nxk))
}

# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
DuchonS <- function(x, m = 2) {
    if (m < 1) 
        stop("'m' has to be larger than zero.")
    if (!is.matrix(x)) 
        x <- as.matrix(x)
    d <- ncol(x)
    n <- nrow(x)
    nterms<- choose((m + d -1),d)
    temp <- .C("polynom", m = as.integer(m), n = as.integer(n), 
        dim = as.integer(d), des = as.double(x), lddes = as.integer(n), 
        npoly = as.integer(nterms), tmatrix = double(n * nterms),
               ldt = as.integer(n), wptr = integer(d * m),
               info = as.integer(0), ptab = integer(nterms * d),
               ldptab = as.integer(nterms))
    temp2 <- matrix(temp$tmatrix, nrow = n)
    return(temp2)
}
