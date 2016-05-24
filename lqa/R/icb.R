icb <-
function (lambda = NULL, ...)
{
   lambda.check (lambda)

   if ((length (lambda) > 2) | (length (lambda) < 2))
     stop ("The genet penalty must consist on two parameters! \n")


   lambda1 <- lambda[1]
   lambda2 <- lambda[2]

   names (lambda) <- c ("lambda1", "lambda2")


  getpenmat <- function (beta = NULL, x = NULL, c1 = lqa.control()$c1, eps.tol = 1e-07, eps.tol2 = 1e-16, ...)
  {
    if (is.null (beta))
      stop ("'beta' must be the current coefficient vector \n")

    l <- list (...)

    if (is.null (x))
      x <- get ("x", envir = l$env)

    x <- as.matrix (x)
    p <- ncol (x)

    penmat <- cor.xv <- matrix (0, nrow = p, ncol = p)
    var.x <- apply (x, 2, var)
    v <- which (var.x != 0)
  
    if (length (v) > 1)
    {
      nv <- which (var.x == 0)
      cor.xv[v,v] <- cor (x[,v])
      cor.xv[v,v] <- cor.xv[v,v] - eps.tol2 * (cor.xv[v,v] == 1) + eps.tol2 * (cor.xv[v,v] == -1)
      penmat[v,v] <- -1 * cor.xv[v,v] / (1 - cor.xv[v,v]^2)
     
      for (i in v)
      {
        penmat[i,i] <- sum (1 / (1 - cor.xv[-c (i, nv), i]^2))
      } 
    }

    if (length (v) == 1)
      penmat[v,v] <- 1

    2 * lambda2 * penmat + lambda1 * diag (1 / (sqrt (beta^2 + c1))) * as.integer (beta != 0)
  }

  structure (list (penalty = "icb", lambda = lambda, getpenmat = getpenmat), class = "penalty")

}

