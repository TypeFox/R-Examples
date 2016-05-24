# extract outer weights from sempls object
plsWeights <- function(object){
  UseMethod("plsWeights", object)
}

plsWeights.sempls <- function(object)
{
  plsWeights <- object$outer_weights
  class(plsWeights) <- c("plsWeights")
  return(plsWeights)
}

print.plsWeights <- function(x, na.print=".", digits=2, abbreviate=FALSE, ...)
{
  weights <- x
  weights[weights==0] <- NA
  if(abbreviate) dimnames(weights) <- lapply(dimnames(weights), abbreviate, ...)
  print.table(weights, na.print=na.print, digits=digits, ...)
  invisible(x)
}
