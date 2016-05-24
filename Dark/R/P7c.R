P7c <- function(a, X) {
  # Seven Param Exp + Bilinear and two fixed transition rates
  x <- .GlobalEnv$x
  y <- .GlobalEnv$y # explicitly locating the values
  if (missing(X)) {
    if (length(a) == 1) {
      list(Pn = 7L, Mod = "P7c")
    } else {
      Yest <- a[1] + a[2] * exp(-x * a[3]) + a[4] * (x - a[5]) * H(x, 10, a[5]) + a[6] * (x - a[7]) * H(x, 
        10, a[7])
      sum((y - Yest) ^ 2)
    }
  } else {
    a[1] + a[2] * exp(-X * a[3]) + a[4] * (X - a[5]) * H(X, 10, a[5]) + a[6] * (X - a[7]) * H(X, 10, a[7])
  }
}
