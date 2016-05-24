# 1 Sample truncated univariate normal vectorized
rtnorm <- function (mu = 0, sd = 1, a = -Inf, b = Inf)
{
  F <- runif(n=length(mu))
  Fa <- pnorm((a - mu)/sd, 0, sd = 1)
  Fa[a == -Inf] <- 0
  Fb <- pnorm((b - mu)/sd, 0, sd = 1)
  Fb[b == Inf] <- 1
  y <- mu + sd * qnorm(F * (Fb - Fa) + Fa)
  y
}