derivchk <- function (x, y, Dy) {

  #  checks that DY is the derivative of Y by comparing it
  #  with Y's central difference estimate.
  #  The value of |DYHAT-DY|/|DY| is returned.

  n <- length(x)
  if (n < 3) stop("X does not have enough elements")
  if (n != length(y) | n != length(Dy)) stop(
        "Lengh of Y or DY not consistent with length of X")
  indup <- 3:n
  inddn <- 1:(n-2)
  indct <- 2:(n-1)
  xdiff <- x[indup]-x[inddn]
  if (min(xdiff) <= 0) stop("X not strictly increasing")
  Dyhat <- (y[indup]-y[inddn])/xdiff
  ratio <- sqrt(mean((Dyhat-Dy[indct])^2))/sqrt(mean(Dy[indct]^2))
  return(ratio)
}
