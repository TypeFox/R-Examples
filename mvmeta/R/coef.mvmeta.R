###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
coef.mvmeta <- 
function(object, format=c("vector","matrix"), ...) {
#
################################################################################
#
  coef <- object$coefficients
  format <- match.arg(format,c("vector","matrix"))
  if(format=="matrix" || is.vector(coef)) return(coef)
  names <- paste(rep(colnames(coef),each=nrow(coef)),
    rep(rownames(coef),ncol(coef)),sep=".")
  coef <- as.numeric(coef)
  names(coef) <- names
#
  coef
}
