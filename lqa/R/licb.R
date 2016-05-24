licb <-
function (lambda = NULL, ...)
{
  lambda.check (lambda)  

  if (length (lambda) != 2)
    stop ("The LICB penalty must consist on two parameters! \n")

  names (lambda) <- c ("lambda1", "lambda2")


  first.derivative <- function (beta = NULL, x = x, ...)   
  {
    if (is.null (beta))
      stop ("'beta' must be the current coefficient vector \n")

    p <- length (beta)
 
    if (p < 2)
      stop ("There must be at least two regressors! \n")

    vec1 <- c (rep (lambda[1], p), rep (lambda[2], (p * (p - 1))))
# Achtung: der 2. Teil enthaelt tatsaechlich p * (p-1) viele Terme, weil wir die Summe aufsplitten!
    return (vec1 * as.integer (drop (t (a.coefs (beta = beta, x = x, ...)) %*% beta) != 0))
  }


  a.coefs <- function (x = NULL, eps.tol2 = 1e-16, ...)
  {
    l <- list (...)

    if (is.null (x))
      x <- get ("x", envir = l$env)

    x <- as.matrix (x)
    p <- ncol (x)

    if (p < 2)
      stop ("There must be at least two regressors! \n")

    var.x <- apply (x, 2, var)
    v <- which (var.x != 0)
    Rmat <- cor (x)
    Rmat[which (Rmat == 1, arr.ind = TRUE)] <- 1 + eps.tol2
    Rmat[which (Rmat == -1, arr.ind = TRUE)] <- -1 + eps.tol2
   
    seq1 <- 1 : (p-1)
    rho.vec <- sapply (seq1, function (seq1) {Rmat[seq1, (seq1 + 1) : p]})

    Amat <- Bmat <- matrix (0, nrow = p, ncol = 0.5 * p * (p - 1))
    l1 <- 0

    for (k in 1 : (p - 1))
    {
       a.vec <- 1 / (1 - rho.vec[[k]])
       b.vec <- 1 / (1 + rho.vec[[k]])
       Amat[k : p, (l1 + 1) : (l1 + p - k)] <- rbind (a.vec, diag (-a.vec, nrow = length (a.vec)))
       Bmat[k : p, (l1 + 1) : (l1 + p - k)] <- rbind (b.vec, diag (-b.vec, nrow = length (b.vec)))
       l1 <- l1 + p - k
    }

    a.coefs.mat <- cbind (diag (p), Amat, Bmat)

    return (a.coefs.mat)
  }


  structure (list (penalty = "licb", lambda = lambda, first.derivative = first.derivative, a.coefs = a.coefs), class = "penalty")
}

