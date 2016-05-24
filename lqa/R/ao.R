ao <-
function (lambda = NULL, ...)
{
   lambda.check (lambda)

   if ((length (lambda) > 2) | (length (lambda) < 2))
     stop ("The ao penalty must consist on two parameters! \n")

   lambda1 <- lambda[1]
   gamma <- lambda[2]

   if (gamma <= 1)
     stop ("ao penalty: 'gamma' must be >= 1! \n")

   names (lambda) <- c ("lambda", "gamma")


  getpenmat <- function (beta = NULL, x = NULL, c1 = lqa.control()$c1, eps.tol = 1e-07, eps.tol2 = 1e-16, ...)
  {
    if (is.null (beta))
      stop ("'beta' must be the current coefficient vector \n")

    l <- list (...)

    if (is.null (x))
      x <- get ("x", envir = l$env)

    x <- as.matrix (x)
    p <- ncol (x)

    if (p < 2)
      stop ("There must be at least two regressors! \n")

    var.x <- apply (x, 2, var)
    v <- which (var.x != 0)
    Rmat <- abs (cor (x))
    seq1 <- 1 : p
    a.vec <- sapply (seq1, function (seq1) {sum (1 - Rmat[seq1, 1 : p])})
    b.vec <- sapply (seq1, function (seq1) {sum (Rmat[seq1, 1 : p]) - 1})

    pen.vec <- (lambda1 * a.vec + gamma * b.vec * abs (beta)^(gamma - 1)) / sqrt (beta^2 + c1) * as.integer (beta != 0)
    return (diag (pen.vec, nrow = p))
  }

  structure (list (penalty = "ao", lambda = lambda, getpenmat = getpenmat), class = "penalty")
}

