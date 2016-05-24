##########################################################################
## function that performs procrustes transformation on X with target Xstar
##
## returns the rotation matrix R, translation vector tt,
##   and dilation factor s for which:
##
##            s X R + 1 tt' \approx Xstar
##
##   along with X.new = s X R + 1 tt'
##
## Based on Borg and Groenen 1997. Modern Multidimensional
##   Scaling. New York: Springer. pp. 340-342.
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Kevin Quinn
## Harvard University
## 6/13/2005
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

procrustes <- function(X, Xstar, translation=FALSE, dilation=FALSE){
  if (nrow(X) != nrow(Xstar)){
    cat("X and Xstar do not have same number of rows.\n")
    stop("Check data and call procrustes() again. \n") 
  }
  if (ncol(X) != ncol(Xstar)){
    cat("X and Xstar do not have same number of columns.\n")
    stop("Check data and call procrustes() again. \n") 
  }

  n <- nrow(X)
  m <- ncol(X)
  
  if (translation){
    J <- diag(n) - 1/n * matrix(1, n, n)
  }
  else{
    J <- diag(n)
  }

  C <- t(Xstar) %*% J %*% X
  svd.out <- svd(C)
  R <- svd.out$v %*% t(svd.out$u)
  s <- 1
  if (dilation){
    mat1 <- t(Xstar) %*% J %*% X %*% R
    mat2 <- t(X) %*% J %*% X
    s.numer <- 0
    s.denom <- 0
    for (i in 1:m){
      s.numer <- s.numer + mat1[i,i]
      s.denom <- s.denom + mat2[i,i]
    }
    s <- s.numer / s.denom
  }
  tt <- matrix(0, m, 1)
  if (translation){
    tt <- 1/n * t(Xstar - s*X %*% R) %*% matrix(1, n, 1)
  }

  X.new <- s * X %*% R + matrix(tt, nrow(X), ncol(X), byrow=TRUE)
  
  return(list(X.new=X.new, R=R, tt=tt, s=s))  
}




