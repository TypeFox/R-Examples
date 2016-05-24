scad <-
function (lambda = NULL, ...)
{
  lambda.check (lambda)

  if (length (lambda) > 2)
    stop ("The scad penalty must consist on two parameters! \n")

  if (length (lambda) == 1)
    lambda[2] <- 3.7    # as suggested in Fan & Li (2001)

  if (lambda[2] <= 2)
    stop ("'lambda[2]' must be '> 2'")

  names (lambda) <- c ("lambda", "a")


  first.derivative <- function (beta = NULL, ...)
  {
    if (is.null (beta))
      stop ("'beta' must be the current coefficient vector \n")

    lambda1 <- lambda[1]
    a <- lambda[2]
    theta <- abs (beta)       # als Argument werden hier die Absolutbeträgt von beta verwendet!!!
     
    p <- length (beta)
    help1 <- sapply (theta, function (theta) {max (a * lambda1 - theta, 0) / ((a - 1) * lambda1)})

    lambda1 * ((theta <= lambda1) + help1 * (theta > lambda1))
  }


  getpenmat <- function (beta = NULL, c1 = lqa.control()$c1, ...)
  {
      if (is.null (beta))
        stop ("'beta' must be the current coefficient vector \n")

      if (c1 < 0)
        stop ("'c1' must be non-negative \n")

      penmat <- diag (first.derivative (beta = beta, ...) / (sqrt (beta^2 + c1) + 0.000001), length (beta))
      penmat
  }


  structure (list (penalty = "scad", lambda = lambda, first.derivative = first.derivative, getpenmat = getpenmat), class = "penalty")
}

