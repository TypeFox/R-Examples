oscar <-
function (lambda = NULL, ...)
{
  lambda.check (lambda)

  if (length (lambda) != 2)
    stop ("The oscar penalty must consist on two parameters! \n")

  names (lambda) = c ("lambda", "c")

  first.derivative <- function (beta = NULL, ...)
  {
    if (is.null (beta))
      stop ("'beta' must be the current coefficient vector \n")

    p <- length (beta)
 
    if (p < 2)
      stop ("There must be at least two regressors! \n")
 
    rep (lambda[1], p) * as.integer (beta != 0)
  }

  a.coefs <- function (beta = NULL, ...)
  {
    if (is.null (beta))
      stop ("'beta' must be the current coefficient vector \n")

    p <- length (beta)
    beta <- abs (beta)
    c2 <- lambda[2]

    order.beta <- order (beta)
    seq1 <- c (c2 * (1 : p - 1) + 1)

    return (diag (seq1[order.beta], length (beta)))
  }

  structure (list (penalty = "oscar", lambda = lambda, first.derivative = first.derivative, a.coefs = a.coefs), class = "penalty")
}

