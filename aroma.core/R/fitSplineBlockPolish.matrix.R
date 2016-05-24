setMethodS3("fitSplineBlockPolish", "matrix", function(z, blockSizes=c(20,20), spar=0.7, maxIter=1, ...) {
  effectFcn <- function(z, x, ...) {
    ok <- (is.finite(x) & is.finite(z));
    xok <- x[ok];
    zok <- z[ok];
    # Not needed anymore
    ok <- NULL;
    fit <- smooth.spline(x=xok, y=zok, ...);
    predict(fit, x=x)$y;
  }

  matrixBlockPolish(z, blockSizes=blockSizes, FUN=effectFcn,
                                           spar=spar, maxIter=maxIter, ...);
})


############################################################################
# HISTORY:
# 2008-04-03
# o BUG FIX: Used NA-filtered 'x' values in predict().
# 2008-04-02
# o Created.  Used to do spatial normalization microarray data.
############################################################################
