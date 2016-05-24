getSigmaR <- function(model, digits=NULL)
{
  ## 1  Parse args
  if(class(model) != "scape")
    stop("The 'model' argument should be a scape object, not ", class(model))

  ## 2  Create output
  output <- model$Dev$sigmaR
  if(!is.null(digits))
    output <- round(output, digits=digits)

  return(output)
}
