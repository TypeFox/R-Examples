bridge <-
function (lambda = NULL, ...)
{
  lambda.check (lambda)  

  if (length (lambda) != 2)
    stop ("The bridge penalty must consist on two parameters! \n")

  if (lambda[2] < 1)
    stop ("'lambda[2]' (= 'gamma') must be '>= 1'")

  lambda1 <- lambda[1]
  gamma <- lambda[2]

  names (lambda) <- c ("lambda", "gamma")


  getpenmat <- function (beta = NULL, c1 = lqa.control()$c1, ...)
  {
      if (is.null (beta))
        stop ("'beta' must be the current coefficient vector \n")

      if (c1 < 0)
        stop ("'c1' must be non-negative \n")

      penmat <- lambda1 * diag (gamma * abs (beta)^(gamma - 1) / (sqrt (beta^2 + c1)), length (beta))
      penmat
  }


  first.derivative <- function (beta, ...)
  {
    p <- length (beta) 
    return (lambda1 * gamma * abs (beta)^(gamma - 1))
  }

  structure (list (penalty = "bridge", lambda = lambda, getpenmat = getpenmat, first.derivative = first.derivative), class = "penalty")
}

