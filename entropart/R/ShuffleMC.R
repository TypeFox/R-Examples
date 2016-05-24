ShuffleMC <- 
function(MClist, Weights = rep(1, length(MClist)), CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()

  # Metacommunities must have names
  if (is.null(names(MClist)))
    names(MClist) <- paste("MC", 1:length(MClist), sep="")
  
  # Merge metacommunities Nsi
  Reduce(function(...) mergeandlabel(...), lapply(MClist, function(x) x$Nsi)) -> Gabundances
  # Suffle the colum order
  Shuffled <- sample(1:ncol(Gabundances))
  # Name the communities
  NumCommunities <- unlist(lapply(MClist, function(x) length(x$Ni)))
  MCnames <- rep(names(MClist), NumCommunities)
  names(Gabundances) <-  paste(MCnames, unlist(lapply(MClist, function(x) names(x$Ni))), sep=".")[Shuffled]
  # Create shuffled MC's
  FirstC <- 1
  ShuffledMCList <- list()
  for (i in 1:length(MClist)) {
    LastC <- FirstC + length((MClist[[i]])$Wi) - 1
    ShuffledMCList[[i]] <- MetaCommunity(Gabundances[, Shuffled[FirstC:LastC]], (MClist[[i]])$Wi)
    names(ShuffledMCList)[i] <- paste("MC", i, sep="")
    FirstC <- LastC + 1
  }
  return(ShuffledMCList)
}
