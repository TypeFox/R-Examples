

rvunif <- function (n=1, min=0, max=1)
{
  rvvapply(stats:::runif, n.=n, min=min, max=max)
}

