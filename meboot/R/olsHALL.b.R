olsHALL.b <- function(y, x)
{
  n <- length(x)
  x <- cbind(1, y[-n], x[-n])
  y <- y[-1]  # matrix(y[-1])
  p <- ncol(x)
  ny <- NCOL(y)
  tol <- 1e-07
  n <- nrow(x)

  #z <- .Fortran("dqrls", qr = x, n = n, p = p, y = y, ny = ny,
  #        tol = as.double(tol), coefficients = mat.or.vec(p, ny),
  #        residuals = y, effects = y, rank = integer(1), pivot = 1:p,
  #        qraux = double(p), work = double(2 * p), PACKAGE = "base")
  z <- lm.fit(x, y)

  z$coefficients
}
