estSigmaR <- function(model, digits=2)
{
  ## 1  Parse args
  if(class(model) != "scape")
    stop("The 'model' argument should be a scape object, not ", class(model))
  if(is.null(model$Dev))
    stop("Element 'Dev' not found; if rec devs were estimated, make sure they were correctly imported")

  ## 2  Calculate sigmaR
  ss <- sapply(model$Dev[c("Initial","Annual")], function(x) sum(x^2))
  n <- sapply(model$Dev[c("Initial","Annual")], length)
  output <- sqrt(ss/(n-1))
  if(!is.null(digits))
    output <- round(output, digits=digits)

  return(output)
}
