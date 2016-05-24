AlphaEntropy <-
function(MC, q = 1, Correction = "Best", Tree = NULL, Normalize = TRUE, Z = NULL, CheckArguments = TRUE) 
{
  if (CheckArguments) 
    CheckentropartArguments()
  
  # Communities
  if (!is.null(Tree)) {
    Method <- "HCDT"
    # Preprocess the tree before apply
    ppTree <- Preprocess.Tree(Tree)
    # Get the entropy of communities.
    DetailedCommunities <- apply(MC$Nsi, 2, bcPhyloEntropy, q=q, Tree=ppTree, Normalize=Normalize, Correction=Correction, CheckArguments=FALSE)
    # Get $Total in each community's list
    Communities <- unlist(lapply(DetailedCommunities, function(x) x$Total))
  } else {
    if (!is.null(Z)) {
      Method <- "Similarity-based"
      Communities <- apply(MC$Nsi, 2, bcHqz, q=q, Correction=Correction, Z=Z, CheckArguments=FALSE)
    } else {
      Method <- "Neutral"
      Communities <- apply(MC$Nsi, 2, bcTsallis, q=q, Correction=Correction, CheckArguments=FALSE)
    }
  }
  # Weighted sum of Communities
  Total <- sum(Communities*MC$Wi)
  
  Entropy <- list(
    MetaCommunity = deparse(substitute(MC)),
    Method = Method,
    Type = "alpha",
    Order = q,
    Correction = Correction,
    Normalized = Normalize,
    Weights = MC$Wi, 
    Communities = Communities, 
    Total=Total
  )
  if(!is.null(Tree))
    Entropy$Tree <- deparse(substitute(Tree)) 
  if(!is.null(Z))
    Entropy$Z <- deparse(substitute(Z)) 
  class(Entropy) <- "MCentropy"
  
  return(Entropy)  
}
