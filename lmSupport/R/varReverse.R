varReverse <- function(Var, LowAnchor, HighAnchor)
{
  if (min(Var, na.rm = TRUE) < LowAnchor) stop('Observed values for variable < LowAnchor')
  if (max(Var, na.rm = TRUE) > HighAnchor) stop('Observed values for variable > HighAnchor')
  
  NewVar = (LowAnchor + HighAnchor) - Var
  
  return(NewVar)
}