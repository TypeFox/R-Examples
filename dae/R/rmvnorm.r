"rmvnorm" <- function(mean, V, method = 'eigenanalysis')
#function to generate a vector of values from an n-dimensional normal 
#distribution whose expectation is given by the n-vector mean and variance by 
#the n x n symmetric matrix V.
#The method:
# a) uses either the eigenvalue or choleski decomposition of the variance matrix 
#    to form a matrix that transforms an iid vector of values to a vector with variance V.
# b) generates a vector of length n of standard normal values
# c) premultiplies the vector of standard normal values by the transpose of the
#    upper triangular factor and, to the result, adds mean.
{ method.opt <- c('eigenanalysis','choleski')
  kopt <- pmatch(method, method.opt)
  if (is.na(kopt))
     stop("method must be either eigenanalysis or choleski")
  n <- length(mean)
  if (!all(dim(V) == c(n, n)))
    stop("Length of mean must equal the number of rows and the number of columns in V")
  if (!isSymmetric(V))
    stop("V must be a symmetric matrix")
  if (method.opt[kopt] == 'eigenanalysis')
  { daeTolerance <- get("daeTolerance", envir=daeEnv)
    #use eigenvalue decomposition to establish transformation matrix
    eigdecomp <- eigen(V, symmetric = TRUE)
    eigenval <- eigdecomp$values
    if (!all(eigenval >= -daeTolerance * abs(max(eigenval))))
       stop("Variance matrix is not nonnegative definite")
    R <- eigdecomp$vectors %*% diag(sqrt(pmax(eigenval, 0)), n)
  }
  else
  { R <- chol(V)
  }
  y <- as.vector(mean + R %*% rnorm(n))
  return(y)
}