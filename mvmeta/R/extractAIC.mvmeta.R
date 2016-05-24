###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
extractAIC.mvmeta <-
function (object, ...) {
#
################################################################################
# EXTRACTS THE NUMBER OF OBSERVATIONS USED FOR FITTING. USED BY BIC
#
   c(object$df$df,AIC(object))
#
}