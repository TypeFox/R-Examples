gaussKernel <- function(x, mean, sd=100e3, wMin=0.3) {
  w <- dnorm(x, mean=mean, sd=sd);
  w <- w / max(w);
  w[w < wMin] <- 0;
  w;
} # gaussKernel()


############################################################################
# HISTORY:
# 2007-09+18
# o Created from gaussianSmoothing.R.
############################################################################
