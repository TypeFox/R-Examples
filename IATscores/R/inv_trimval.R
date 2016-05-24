inv_trimval <- function (x, tr = .1)
{
  #inv_trimval is a function that trims central values
  # it does NOT accept vectors with NAs, use trim_or_win()
  y <- sort(x)
  n <- length(x)
  ibot <- ceiling((.5 - tr) * n) + 1
  itop <- n - ibot + 1
  xbot <- y[ibot]
  xtop <- y[itop]
  x[x >= xbot & x <= xtop] <- NA
  x
}