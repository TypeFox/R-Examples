varRecode <- function(Var, Old, New)
{
  #check for length match for levels of old and new
  if (length(Old) != length(New))
    stop('Number of Old and New values do not match')

  NewVar = rep(NA, length(Var)) #set to NA to start
  for (i in 1:length(Old))
  {
    NewVar[Var==Old[i]] = New[i]
  }
  
  #check for mismatch of NA which means some levels in Var were not reassigned.  Warning only
  if (sum(is.na(NewVar)) != sum(is.na(Var)))
    warning(sprintf('%.0f NAs in original variable but %.0f NAs in recoded variable', sum(is.na(Var)), sum(is.na(NewVar))))

  #covert NewVar to factor if Var was factor
  if (is.factor(Var))  
    NewVar = as.factor(NewVar)
  
  return(NewVar)
}