rRandomLabeling <-
function(X, CheckArguments = TRUE) {
  
  if (CheckArguments)
    CheckdbmssArguments()
  
  return (rRandomLocation(X, CheckArguments=FALSE))
}