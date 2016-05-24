
dgp <- function(n, sd = 1) {
   x1 <- seq(from = 0, to = 2 * pi, length = n)
   x2 <- runif(n = n)
   y <- rnorm(n = n, mean = sin(x1) + x2, sd = sd)
   return(data.frame(y, x1, x2))
}
