#' @export
getCall <- function(MLCall){
  # Get the functionname
  bracketLocation <- gregexpr(pattern ='\\(',MLCall)[[1]][1]
  parentFunction  <- substr(MLCall,1,bracketLocation-1)
  
  # Parse the call to a function call
  functionCall    <- match.call(eval(parse(text=parentFunction)),parse(text=MLCall))
  
  # Return the function call
  return(functionCall)
}