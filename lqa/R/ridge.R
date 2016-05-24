ridge <-
function (lambda = NULL, ...)
{
  lambda.check (lambda)  

  if (length (lambda) != 1)        ## Check on dimensionality of lambda
    stop ("lambda must be a scalar \n")

  names (lambda) <- "lambda"

  getpenmat <- function (beta = NULL, ...)
  {
      if (is.null (beta))
        stop ("'beta' must be the current coefficient vector \n")

      p <- length (beta)
      penmat <- diag (lambda, p)
      penmat
  }

  structure (list (penalty = "ridge", lambda = lambda, getpenmat = getpenmat), class = "penalty")
}

