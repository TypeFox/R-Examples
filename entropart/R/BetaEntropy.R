BetaEntropy <-
function(MC, q = 1, Correction = "Best", Tree = NULL, Normalize = TRUE, Z = NULL, CheckArguments = TRUE) 
{
  if (CheckArguments) 
    CheckentropartArguments()
    
  # Communities
  if (!is.null(Tree)) {
    Method <- "HCDT"
    # Get the entropy of communities. Unlist to get a vector
    DetailedCommunities <- apply(MC$Nsi, 2, bcPhyloBetaEntropy, Nexp=MC$Ns, q=q, Tree=Tree, Normalize=Normalize, Correction=Correction, CheckArguments=FALSE)
    # Get $Total in each community's list
    Communities <- unlist(lapply(DetailedCommunities, function(x) x$Total))
  } else {
    if (!is.null(Z)) {
      Method <- "Similarity-based"
      Communities <- apply(MC$Nsi, 2, bcHqzBeta, Nexp=MC$Ns, q=q, Correction=Correction, Z=Z, CheckArguments=FALSE)
    } else {
      Method <- "Neutral"
      Communities <- apply(MC$Nsi, 2, bcTsallisBeta, Nexp=MC$Ns, q=q, Correction=Correction, CheckArguments=FALSE)
    }
  }
  # Weighted sum of Communities
  Total <- sum(Communities*MC$Wi)
  
  Entropy <- list(
    MetaCommunity = deparse(substitute(MC)),
    Method = Method,
    Type = "beta",
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
