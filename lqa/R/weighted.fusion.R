weighted.fusion <-
function (lambda = NULL, ...)
{
   lambda.check (lambda)

   if ((length (lambda) > 3) | (length (lambda) < 2))
     stop ("The weighted fusion penalty must consist on two (enet) or three parameters (genet)! \n")

   if (length (lambda) == 2)
     lambda[3] <- 2

   lambda1 <- lambda[1]
   lambda2 <- lambda[2]
   gamma <- lambda[3]

   if (gamma <= 0)
     stop ("'lambda[3]' (= gamma) must be greater than zero \n")

   names (lambda) <- c ("lambda1", "lambda2", "gamma")


  getpenmat <- function (beta = NULL, c1 = lqa.control()$c1, x, cor.mat, ...)
  {
      if (is.null (beta))
        stop ("'beta' must be the current coefficient vector \n")

      if (c1 < 0)
        stop ("'c1' must be non-negative \n")

      if (missing (x))
      {
        l <- list (...)
        x <- get ("x", envir = l$env)
      }

      if (var (x[,1]) <= lqa.control()$var.eps)   # check (and correct) for intercept
        x <- x[,-1]

      if (missing (cor.mat))
        cor.mat <- cor (x)

      lasso.part <- lambda1 * diag (1 / sqrt (beta^2 + c1))
     
      abs.cor.mat <- abs (cor.mat)
      sign.mat <- sign (cor.mat) - diag (2, nrow (cor.mat))  # includes minus ones at the diagonal
      Wmat <- (abs.cor.mat^gamma - 1 * (abs.cor.mat == 1)) / (1 - abs.cor.mat * (abs.cor.mat != 1))   # this sets the diagonal elements to zero
      Wj.vec <- apply (Wmat, 1, sum)
      fusion.part <- lambda2 * -sign.mat * (Wmat + diag (Wj.vec))
      
      penmat <- lasso.part + fusion.part
      penmat
  }

  structure (list (penalty = "weighted.fusion", lambda = lambda, getpenmat = getpenmat), class = "penalty")
}

