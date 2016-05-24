genet <-
function (lambda = NULL, ...)
{
   lambda.check (lambda)

   if ((length (lambda) > 3) | (length (lambda) < 2))
     stop ("The genet penalty must consist on two (enet) or three parameters (genet)! \n")

   if (length (lambda) == 2)
     lambda[3] <- 2

   lambda1 <- lambda[1]
   alpha <- lambda[2]
   gamma <- lambda[3]

   if (alpha > 1) 
     stop ("'lambda[2]' (= alpha) must be between 0 and 1\n")

   if (gamma <= 1)
     stop ("'lambda[3]' (= gamma) must be greater than one \n")

   names (lambda) <- c ("lambda", "alpha", "gamma")


  getpenmat <- function (beta = NULL, c1 = lqa.control()$c1, ...)
  {
      if (is.null (beta))
        stop ("'beta' must be the current coefficient vector \n")

      if (c1 < 0)
        stop ("'c1' must be non-negative \n")

      penmat <- lambda1 * diag (((1 - alpha) * gamma * abs (beta)^(gamma - 1) + alpha) / (sqrt (beta^2 + c1)), length (beta))
      penmat
  }

  structure (list (penalty = "genet", lambda = lambda, getpenmat = getpenmat), class = "penalty")
}

