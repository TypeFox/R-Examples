adaptive.lasso <-
function (lambda = NULL, al.weights = NULL, ...)
{
  lambda.check (lambda)  

  if (length (lambda) != 1)        ## Check on dimensionality of lambda
    stop ("lambda must be a scalar \n")

  names (lambda) <- "lambda"

  if (is.null (al.weights))
    al.weights <- 1


  getpenmat <- function (beta = NULL, c1 = lqa.control()$c1, ...)
  {
      if (is.null (beta))
        stop ("'beta' must be the current coefficient vector \n")
  
      if (c1 < 0)
        stop ("'c1' must be non-negative \n")

      penmat <- lambda * diag (al.weights / (sqrt (beta^2 + c1))) * as.integer (beta != 0)
      penmat
  }


  first.derivative <- function (beta, ...)
  {
    p <- length (beta) 
    return (rep (lambda * al.weights, p))
  }

  structure (list (penalty = "adaptive.lasso", lambda = lambda, getpenmat = getpenmat, first.derivative = first.derivative), class = "penalty")
}

