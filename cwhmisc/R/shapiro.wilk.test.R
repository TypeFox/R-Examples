"shapiro.wilk.test" <- function(x){
  ##this function is an S version of the procedure described by
  # J. P. Royston (1982) in "An Extension of Shapiro and Wilk's W Test
  # for Normality to Large Samples" from Applied Statistics,31 no.2
  #  pp115-124.
  #
  n   <- length(x)
  index <- 1:n
  m   <- qnorm((index - 0.375)/(n + 0.25))
  y   <- sort(x)
  mu  <- mean(y)
  SSq <- sum((y - mu)^2)
  astar <- 2 * m
  ends  <- c(1, n)
  astar.p <- astar[ - ends]
  if (n <= 20) m <- n - 1
  else m <- n
  if (m < 20) aa <- gamma(0.5 * (m + 1))/(sqrt(2) * gamma(0.5 * m + 1))
  else {
    f1 <- (6 * m + 7)/(6 * m + 13)
    f2 <- exp(1)/(m + 2)
    f3 <- (m + 1)/(m + 2)
    f3 <- f3^(m - 2)
    aa <- f1 * sqrt(f2 * f3)
  }
  astar.1  <- (aa * sum(astar.p^2))/(1 - 2 * aa)
  astar.1  <- sqrt(astar.1)
  astar[1] <-  - astar.1
  astar[n] <- astar.1
  A <- astar/sqrt(sum(astar^2))
  W <- (sum(A * y)^2)/SSq
  if (n <= 20) {
    u <- log(n) - 3
    lambda <- 0.118898 + 0.133414 * u + 0.327907 * u^2
    logmu <- -0.37542 - 0.492145 * u - 1.124332 * u^2 - 0.199422 * u^3
    logsigma <- -3.15805 + 0.729399 * u + 3.01855 * u^2 + 1.558776 * u^3
  } else {
    u <- log(n) - 5
    lambda <- 0.480385 + 0.318828 * u - 0.0241665 * u^3 + 0.00879701 * u^4 + 0.002989646 * u^5
    logmu <- -1.91487 - 1.37888 * u - 0.04183209 * u^2 + 0.1066339 * u^3 - 0.03513666 * u^4 - 0.01504614 * u^5
    logsigma <- -3.73538 - 1.015807 * u - 0.331885 * u^ 2 + 0.1773538 * u^3 - 0.01638782 * u^4 - 0.03215018 * u^5 + 0.003852646 * u^6
  }
  mu <- exp(logmu)
  sigma <- exp(logsigma)
  y <- (1 - W)^lambda
  z <- (y - mu)/sigma
  p <- 1 - pnorm(z)
  if(n < 7) {
    warning("n is too small for this program to correctly estimate p")
    p <- NA
  } else  if (n > 2000) {
    warning("n is too large for this program to correctly estimate p")
    p <- NA
  }
  list(W = W, n = n, p = p)
}
