if (!isGeneric("train")) {
  setGeneric("train")
}
if (!isGeneric("featurePlot")) {
  setGeneric("featurePlot")
}

setMethod("featurePlot", signature(x = ".CaretHyperspectral"),
          definition = function(x,
                                y,
                                ...)
{
  y_missing <- missing(y)
  
  if (y_missing)
  {
    y <- .getResponseVar(x, 
                         advice = c("train", "setResponse", 
                                    "This is only required if you do not specify 'y'.")) 
  }
  
  useAttributesAsPredicants <- !is.na(.getPredicantVar(x, stopifmissing = FALSE))[1]
  
  all_vals <- as.data.frame(x)
  
  if (useAttributesAsPredicants)
  {
    addVar <- .getPredicantVar(x)
    all_vals <- cbind(all_vals, addVar)
  }
  
  x <- all_vals
  return(featurePlot(x = x, y = y, ...))
})
