###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
nobs.mvmeta <-
function (object, ...) {
#
################################################################################
# EXTRACTS THE NUMBER OF OBSERVATIONS USED FOR FITTING. USED BY BIC
#
  object$df$nobs
#
}