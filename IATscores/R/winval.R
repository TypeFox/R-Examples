winval <- function (x, tr = 0.1) 
{
  # copy-pasted from package "WRS"
  # it does NOT accept vectors with NAs, use trim_or_win()
  
  y <- sort(x)
  n <- length(x)
  ibot <- floor(tr * n) + 1
  itop <- n - ibot + 1
  xbot <- y[ibot]
  xtop <- y[itop]
  winval <- ifelse(x <= xbot, xbot, x)
  winval <- ifelse(winval >= xtop, xtop, winval)
  winval
}