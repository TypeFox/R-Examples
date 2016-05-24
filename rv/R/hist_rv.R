# ========================================================================
# hist.rv  -  histogram, adapted for rv's
# ========================================================================

hist.rv <- function(x, grid=c(4,5), xlim=x.range, main=paste(xname,"simulation"), freq=FALSE, ...) {
  par(mfrow=grid)
  #  l <- length(x)
  #  par(mfrow=c(l %% 3 + l %/% 3, 3))
  xname <- deparse(substitute(x))
  grid.n <- grid[1]*grid[2]
  s <- sims(x)
  x.range <- c(min(s),max(s))
  if (grid.n<1)
    stop("Bad grid")
  s <- s[1:grid.n,]
  for (i in 1:grid.n) {
    # truehist(s[i,], xlim=xlim, main, ...)
    hist(s[i,], xlim=xlim, main=main, freq=freq, ...)
  }
}

