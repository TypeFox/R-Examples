setMethodS3("backtransformGenotypeCone", "matrix", function(y, fit, ...) {
  # x <- t(t(fit$Winv) %*% (t(y) - fit$origin));
  x <- t(t(y) - fit$origin) %*% fit$Winv;
  attr(x, "fit") <- fit;
  x;
})


############################################################################
# HISTORY:
# 2006-05-08
# o Created.
############################################################################
