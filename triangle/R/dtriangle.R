################################################################################
#
# Program:   dtriangle.R
# Purpose:   To calculate the PDF for the triangle distribution
# Author:    Rob Carnell
# Date:      June 05
#
# Variables
#   used the same naming conventions as other R distributions (q,p,d)
#   p = cumulative probability
#   a = left triangle endpoint,
#   b = right triangle endpoint
#   c = distribution mode
#   First, exclude situations which are impossible with the function definition
#   Next, define the value of the function on the various intervals
#
################################################################################

dtriangle <- function(x, a=0, b=1, c=(a+b)/2) {
  x1 <- length(x)
  a1 <- length(a)
  b1 <- length(b)
  c1 <- length(c)

  dTest <- function(X){
    if(any(is.na(X))){ # is.na is TRUE for NA, NaN, and FALSE
      if(any(is.nan(X))) return(NaN) # to conform to qunif
      else return(NA) # to conform to qunif
    } else if(X[2] > X[4] | X[3] < X[4] | (X[1]==X[2] & X[2]==X[4])){
      warning("values required to be  a <= c <= b (at least one strict inequality)")
      return(NaN) # to conform to behavior of qunif
    } else if(any(is.infinite(X[2:4]))){
      return(NaN)
    } else if(X[1] <= X[2]) {
      return(0)
    } else if(X[2] != X[4] & X[1] < X[4]){
      return(2*(X[1] - X[2]) / (X[3] - X[2]) / (X[4] - X[2]))
    } else if(X[4] != X[3] & X[1] >= X[4] & X[1] < X[3]){
      return(2*(X[3] - X[1]) / (X[3] - X[2]) / (X[3] - X[4]))
    } else if(X[1] >= X[3]){
      return(0)
    }
  }

  k <- max(x1, a1, b1, c1)
  if(k==1) return(dTest(c(x, a, b, c)))

  params <- matrix(nrow=k, ncol=4)
  tryCatch(
  {
    params[,1] <- x
    params[,2] <- a
    params[,3] <- b
    params[,4] <- c
  }, error = function(X) {
    stop(paste(" -- Argument Lengths: length of x = ", x1,
                ", a = ", a1, ", b = ", b1, ", c = ", c1, " -- ", X, sep=""))
  })

  return(apply(params, 1, dTest))
}

