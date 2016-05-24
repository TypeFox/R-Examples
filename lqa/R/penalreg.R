penalreg <-
function (lambda = NULL, ...)
{
  lambda.check (lambda)
 
  if (length (lambda) != 1)        ## Check on dimensionality of lambda
    stop ("lambda must be a scalar \n")

  names (lambda) <- "lambda"

  getpenmat <- function (x = NULL, eps.tol = 1e-07, eps.tol2 = 1e-16, ...)
  {
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

    4 * lambda * penmat
  }

  structure (list (penalty = "penalreg", lambda = lambda, getpenmat = getpenmat), class = "penalty")

}

