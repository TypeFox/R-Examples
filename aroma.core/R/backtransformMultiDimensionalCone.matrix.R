setMethodS3("backtransformMultiDimensionalCone", "matrix", function(y, fit, ...) {
  # x <- t(t(fit$Winv) %*% (t(y) - fit$origin));
  x <- t(t(y) - fit$origin) %*% fit$Winv;
  attr(x, "fit") <- fit;
  x;
})


############################################################################
# HISTORY:
# 2008-08-02
# o Created from backtransformGenotypeCone.matrix.R.
############################################################################
