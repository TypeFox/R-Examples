GenSeries <- function(n, alpha, mean, std) {
  res <- vector("numeric", n)
  x <- mean
  stdterm <- std * (sqrt(1 - alpha ^ 2) / (1 - alpha))
  for (i in 1:100) {
    x <- alpha * x + (1 - alpha) * rnorm(1, mean, stdterm)
  }
  for (i in 1:n) {
    x <- alpha * x + (1 - alpha) * rnorm(1, mean, stdterm)
    res[i] <- x
  }
  
  res
}
