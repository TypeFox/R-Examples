FitAcfCoef <- function(a, b) {
  if ((a < 0) || (a > 1) || (b < 0) || (b > 1)) {
    stop("One of the coefficients is outside the [0 1] interval");
  }
  p <- (2 - b)
  q <- -2 * a
  delta <- q**2 + 4 / 27 * p**3

  minimum <- (((-q + sqrt(delta)) * 0.5) ^ (1 / 3) 
              - ((q + sqrt(delta)) * 0.5) ^ (1 / 3))
          
  minimum 
}
