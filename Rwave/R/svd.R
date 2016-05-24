
#**************************************************************#
#              (c) Copyright  1997                             #
#                         by                                   #
#     Author: Rene Carmona, Bruno Torresani, Wen L. Hwang      #
#                 Princeton University                         #
#                 All right reserved                           #
#**************************************************************#


## compute the svd (single value decomposition)
## The Input/Output of the command is the same as the commands svd in Splus

SVD <- function(a) {
  m <- dim(a)[1]
  n <- dim(a)[2]
  
  ## diagonal elements
  w <- numeric(n)
  dim(w) <- c(n,1)	 
  
  v <- matrix(0,n,n)
  v <- c(v)
  dim(v) <- c(length(v),1)
  a <- c(a)
  dim(a) <- c(length(a),1)
  
  z <- .C("Ssvdecomp",
          u= as.double(a),
          as.integer(m),
          as.integer(n),
          w = as.double(w),
          v = as.double(v),
          PACKAGE="Rwave")
  
  u <- z$u
  dim(u) <- c(n,n)	
  w <- z$w
  dim(w) <- c(n,1)
  v <- z$v
  dim(v) <- c(n,n)
  
  list(d=w, u=u, v=v)
}





