BetaDiversity <-
function(MC, q = 1, Correction = "Best", Tree = NULL, Normalize = TRUE, Z = NULL, CheckArguments = TRUE) 
{
  if (CheckArguments) 
    CheckentropartArguments()
  
  # Preprocess the tree for performance
  ppTree <- Preprocess.Tree(Tree)
  
  BetaEntropy  <- BetaEntropy (MC=MC, q=q, Correction=Correction, Tree=ppTree, Normalize=Normalize, Z=Z, CheckArguments=FALSE)
  AlphaEntropy <- AlphaEntropy(MC=MC, q=q, Correction=Correction, Tree=ppTree, Normalize=Normalize, Z=Z, CheckArguments=FALSE)

  Diversity <- list(
    MetaCommunity = deparse(substitute(MC)),
    Method = AlphaEntropy$Method,
    Type = "beta",
    Order = q,
    Correction = Correction,
    Normalized = Normalize,
    Weights = MC$Wi, 
    Total = expq(BetaEntropy$Total / (1 - (q-1)*AlphaEntropy$Total), q)
  )
  if(!is.null(Tree))
    Diversity$Tree <- deparse(substitute(Tree)) 
  if(!is.null(Z))
    Diversity$Z <- deparse(substitute(Z)) 
  class(Diversity) <- "MCdiversity"
  
  return(Diversity)  
  
}
