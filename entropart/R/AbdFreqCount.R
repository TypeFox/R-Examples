AbdFreqCount <- 
function (Ns, CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()
  
  NsInt <- as.integer(round(Ns))
  if (any(abs(NsInt-Ns) > sum(Ns)*.Machine$double.eps)) warning ("The abundance frequency count requires integer abundances. Abundances have been rounded.")
  
  # Eliminate 0
  Ns <- NsInt[NsInt > 0]
  
  DistNs <- tapply(Ns, Ns, length)
  afc <- matrix(c(as.integer(names(DistNs)), DistNs), ncol = 2)
  colnames(afc) <- c("Abundance", "NbSpecies")
  return(afc)
}