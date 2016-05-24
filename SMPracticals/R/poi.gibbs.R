"poi.gibbs" <-
function(d, alpha, gamma, delta, I, S)
{
  poi.sim <- function(y, x, alpha, gamma, delta, theta)
  {
    n <- length(y)
    lambda <- theta[1:n]
    beta <- theta[n + 1]
    out1 <- rgamma(n, alpha + y)/(x + 1/beta)
    out2 <- (sum(out1) + delta)/rgamma(1, gamma + n * alpha)
    c(out1, out2)
  }
  n <- length(d$y)
  out <- array(NA, dim = c(I, S,n+1))
  out[1, ,  ] <- rexp((n + 1) * S)
  for(s in 1:S)
  for(i in 2:I)
   out[i,s,] <- poi.sim(d$y, d$x, alpha, gamma, delta, out[i - 1, s,])
  out
}

