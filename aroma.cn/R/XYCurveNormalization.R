###########################################################################/**
# @RdocClass XYCurveNormalization
#
# @title "The XYCurveNormalization class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @see "AbstractCurveNormalization".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("XYCurveNormalization", function(...) {
  extend(AbstractCurveNormalization(...), "XYCurveNormalization");
})

setMethodS3("fitOne", "XYCurveNormalization", function(this, theta, ...) {
  .fitXYCurve(theta, ...);
}, protected=TRUE)

setMethodS3("backtransformOne", "XYCurveNormalization", function(this, theta, fit, ...) {
  .backtransformXYCurve(theta, fit=fit, ...);
}, protected=TRUE)



############################################################################
# HISTORY:
# 2009-07-15
# o Created.
############################################################################
