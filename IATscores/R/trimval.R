trimval <- function (x, tr = .1)
{
  # trimval is a function that trims extreme values
  # it does NOT accept vectors with NAs, use trim_or_win()
  y <- sort(x)
  n <- length(x)
  ibot <- floor(tr * n) + 1
  itop <- n - ibot + 1
  xbot <- y[ibot]
  xtop <- y[itop]
  x[x < xbot] <- NA
  x[x > xtop] <- NA
  x
}