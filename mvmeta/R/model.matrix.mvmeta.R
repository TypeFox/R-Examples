###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
model.matrix.mvmeta <-
function(object, ...) {
#
################################################################################
#
  data <- model.frame(object,xlev=object$xlevels, ...)
  NextMethod("model.matrix",data=data,contrasts.arg=object$contrasts)
}