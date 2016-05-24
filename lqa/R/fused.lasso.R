fused.lasso <-
function (lambda = NULL, ...)
{
  lambda.check (lambda)  

  if (length (lambda) != 2)
    stop ("The fused.lasso penalty must consist on two parameters! \n")

  names (lambda) <- c ("lambda1", "lambda2")


  first.derivative <- function (beta = NULL, ...)   # returns d p_{\lambda,j}(|\a_j^\top \betafat|) / d |\a_j^\top \betafat|  
  {
    if (is.null (beta))
      stop ("'beta' must be the current coefficient vector \n")

    p <- length (beta)
 
    if (p < 2)
      stop ("There must be at least two regressors! \n")

    vec1 <- c (rep (lambda[1], p), rep (lambda[2], p - 1))
    return (vec1 * as.integer (drop (t (a.coefs (beta)) %*% beta) != 0))
  }

  a.coefs <- function (beta = NULL, ...)
  {
    if (is.null (beta))
      stop ("'beta' must be the current coefficient vector \n")

    p <- length (beta) 
 
    if (p < 2)
      stop ("There must be at least two regressors! \n")

    if (p > 2)
    {
      h1 <- cbind (-diag (p-1), 0)
      h2 <- cbind (0, diag (p-1))
      mat1 <- h1 + h2
      mat2 <- diag (p)
      a.coefs.mat <- cbind (mat2, t (mat1))
    }
    else
      a.coefs.mat <- cbind (diag (2), c(-1,1))

    return (a.coefs.mat)
  }


  structure (list (penalty = "fused.lasso", lambda = lambda, first.derivative = first.derivative, a.coefs = a.coefs), class = "penalty")
}

