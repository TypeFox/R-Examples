GammaEntropy <-
function(MC, q = 1, Correction = "Best", Tree = NULL, Normalize = TRUE, Z = NULL, CheckArguments = TRUE) 
{
  if (CheckArguments) 
    CheckentropartArguments()
  
  if (!is.null(Tree)) {
    # Method <- "HCDT"
    Entropy <- bcPhyloEntropy(MC$Ns, q, Tree, Normalize, Correction, CheckArguments=FALSE)$Total
  } else {
    if (!is.null(Z)) {
      # Method <- "Similarity-based"
      Entropy <- bcHqz(MC$Ns, q, Z, Correction, CheckArguments=FALSE)
    } else {
      # Method <- "HCDT"
      Entropy <- bcTsallis(MC$Ns, q, Correction, CheckArguments=FALSE)
    }
  }
  
  return(Entropy)
}
