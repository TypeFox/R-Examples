MergeMC <- 
function(MClist, Weights = rep(1, length(MClist)), CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Metacommunities must have names
  if (is.null(names(MClist)))
    names(MClist) <- paste("MC", 1:length(MClist), sep="")

  # Merge metacommunities Ns
  Reduce(function(...) mergeandlabel(...), lapply(MClist, function(x) x$Ns)) -> Gabundances
  names(Gabundances) <- names(MClist)
  
  # Create the global MC
  return(MetaCommunity(Gabundances, Weights))
}
