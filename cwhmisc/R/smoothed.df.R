smoothed.df <- function(d) {
    F <- cumsum(d$y)
    F <- F/F[length(F)]
    splinefun(d$x, F)
  }
