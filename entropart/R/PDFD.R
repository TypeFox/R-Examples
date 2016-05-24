PDFD <-
function(Ps, Tree, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  return (AllenH(Ps, 0, Tree, Normalize=FALSE, CheckArguments=FALSE))
}
