if (!isGeneric("train")) {
  setGeneric("train")
}
# if (!isGeneric("train.formula")) {
#   setGeneric("train.formula")
# }

setMethod("train", signature(x = ".CaretHyperspectral"),
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
  
  if (class(x) == "Nri")
  {
    all_vals <- as.data.frame(x, na.rm = TRUE)
  } else {
    all_vals <- as.data.frame(x)
  }
  
  if (useAttributesAsPredicants)
  {
    addVar <- .getPredicantVar(x)
    all_vals <- cbind(all_vals, addVar)
  }
  
  x <- all_vals
  return(train(x = x, y = y, ...))
})

# setMethod("train.formula", signature(form = "formula", data = "Speclib"),
#           definition = function(form,
#                                 data,
#                                 ...)
# {
#   if (is.null(bandnames(spectral_data)))
#     bandnames(spectral_data) <- paste("V", c(1:nbands(data)), sep = "")
#   data <- cbind(attribute(data), as.data.frame(data))
#   return(train.formula(form = form, data = data, ...))
# })
