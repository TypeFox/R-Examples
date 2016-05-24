GammaDiversity <-
function(MC, q = 1, Correction = "Best", Tree = NULL, Normalize = TRUE, Z = NULL, CheckArguments = TRUE) 
{
  if (CheckArguments) 
    CheckentropartArguments()
  
  ppTree <- Preprocess.Tree(Tree)
  if (Normalize) {
    Height <- 1
  } else {
    Height <- ppTree$Height
  }  
  GammaEntropy <- GammaEntropy(MC, q, Correction, ppTree, Z, Normalize=TRUE, CheckArguments=FALSE)
  
  return(expq(GammaEntropy, q) * Height)  
}
