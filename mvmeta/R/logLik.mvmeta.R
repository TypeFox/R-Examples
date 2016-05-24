###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
logLik.mvmeta <-
function(object, ...) {
#
################################################################################
#
  val <- object$logLik
  attributes(val) <- object$df
#
  class(val) <- "logLik"
#
  val
}

