hard.rejection <- function(distances, p, beta = 0.9, ...)
{
  d0 <- qchisq(beta, p) * median(distances) / qchisq(0.5, p)
  weights <- double(length(distances))
  weights[distances <= d0] <- 1.0
  weights
}
