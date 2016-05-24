MergeC <- 
function(MClist, Weights = rep(1, length(MClist)), CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Metacommunities must have names
  if (is.null(names(MClist)))
      names(MClist) <- paste("MC", 1:length(MClist), sep="")
  
  CommunityNames <- function(MClist) {
    MCnames <- rep(names(MClist), unlist(lapply(MClist, function(x) length(x$Ni))))
    paste(MCnames, unlist(lapply(MClist, function(x) names(x$Ni))), sep=".")
  }
  # Merge metacommunities Nsi
  Reduce(function(...) mergeandlabel(...), lapply(MClist, function(x) x$Nsi)) -> Gabundances
  NumCommunities <- unlist(lapply(MClist, function(x) length(x$Ni)))
  MCnames <- rep(names(MClist), NumCommunities)
  names(Gabundances) <-  paste(MCnames, unlist(lapply(MClist, function(x) names(x$Ni))), sep=".")
  MCWeights <- unlist(lapply(MClist, function(x) x$Wi))
  
  # Create the global MC
  return(MetaCommunity(Gabundances, MCWeights*rep(Weights, NumCommunities)))
}
