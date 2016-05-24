enet <-
function (lambda = NULL, ...)
{
   lambda.check (lambda)

   if ((length (lambda) > 2) | (length (lambda) < 2))
     stop ("The genet penalty must consist on two parameters! \n")


   lambda1 <- lambda[1]
   lambda2 <- lambda[2]


   names (lambda) <- c ("lambda1", "lambda2")


  getpenmat <- function (beta = NULL, c1 = lqa.control()$c1, ...)
  {
      if (is.null (beta))
        stop ("'beta' must be the current coefficient vector \n")

      if (c1 < 0)
        stop ("'c1' must be non-negative \n")

      penmat <- lambda1 * diag (1 / (sqrt (beta^2 + c1))) * as.integer (beta != 0) + diag (lambda2, length (beta))
      penmat
  }

  structure (list (penalty = "genet", lambda = lambda, getpenmat = getpenmat), class = "penalty")
}

