"BesselPolynomial" <-
function(n, reverse = FALSE)
{
  k <- 0:n
  y.2 <- 1
  y.1 <- c(1, 1)
  if(n == 0) {
    y <- y.2
  } else if(n == 1) {
    y <- y.1
  } else {
    for(i in 2:n) {
      y <- ( 2 * i - 1 ) * c(0, y.1) + c(y.2, 0, 0)
      y.2 <- y.1
      y.1 <- y
    }
  }
  if(reverse) rev(y) else y # if reverse return coefficients from highest to lowest
}
