PhyloEntropy <-
function(NorP, q = 1, Tree, Normalize = TRUE, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  UseMethod("PhyloEntropy")
}


PhyloEntropy.ProbaVector <-
function(NorP, q = 1, Tree, Normalize = TRUE, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ps)) {
      NorP <- Ps
    } else {
      stop("An argument NorP or Ps must be provided.")
    }
  }
  if (CheckArguments)
    CheckentropartArguments()
  
  # Calculate the PhyloValue
  Entropy <- PhyloApply(Tree, Tsallis, NorP, Normalize, q=q, CheckArguments=FALSE)
  # Complete it
  Entropy$Function <- "PhyloEntropy" 
  Entropy$Distribution <- deparse(substitute(NorP)) 
  Entropy$Tree <- deparse(substitute(Tree))
  Entropy$Type <- "alpha or gamma"
  Entropy$Order <- q
  
  class(Entropy) <- c("PhyloEntropy", class(Entropy))
  
  return (Entropy)
}


PhyloEntropy.AbdVector <-
function(NorP, q = 1, Tree, Normalize = TRUE, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return(bcPhyloEntropy(Ns=NorP, q=q, Tree=Tree, Normalize=Normalize, Correction=Correction, CheckArguments=CheckArguments))
}


PhyloEntropy.integer <-
function(NorP, q = 1, Tree, Normalize = TRUE, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return(bcPhyloEntropy(Ns=NorP, q=q, Tree=Tree, Normalize=Normalize, Correction=Correction, CheckArguments=CheckArguments))
}


PhyloEntropy.numeric <-
function(NorP, q = 1, Tree, Normalize = TRUE, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ps)) {
      NorP <- Ps
    } else {
      if (!missing(Ns)) {
        NorP <- Ns
      } else {
        stop("An argument NorP or Ps or Ns must be provided.")
      }
    }
  }
  
  if (abs(sum(NorP) - 1) < length(NorP)*.Machine$double.eps) {
    # Probabilities sum to 1, allowing rounding error
    return(PhyloEntropy.ProbaVector(NorP, q=q, Tree=Tree, Normalize=Normalize, CheckArguments=CheckArguments))
  } else {
    # Abundances
    return(PhyloEntropy.AbdVector(NorP, q=q, Tree=Tree, Normalize=Normalize, Correction=Correction, CheckArguments=CheckArguments))
  }
}


bcPhyloEntropy <-
function(Ns, q = 1, Tree, Normalize = TRUE, Correction = "Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Calculate the PhyloValue
  Entropy <- PhyloApply(Tree, bcTsallis, Ns, Normalize, q=q, Correction=Correction, CheckArguments=FALSE)
  # Complete it
  Entropy$Function <- "PhyloEntropy" 
  Entropy$Distribution <- deparse(substitute(Ns)) 
  Entropy$Tree <- deparse(substitute(Tree))
  Entropy$Type <- "alpha or gamma"
  Entropy$Order <- q
  Entropy$Correction <- Correction
  
  class(Entropy) <- c("PhyloEntropy", class(Entropy))
  
  return (Entropy)                         
}


is.PhyloEntropy <-
function (x) 
{
  inherits(x, "PhyloEntropy")
}


summary.PhyloEntropy <-
function(object, ...) 
{
  cat(object$Type, "phylogenetic or functional entropy of order", object$Order, "of distribution", object$Distribution, fill=TRUE)
  if (!is.null(object$Correction)) {
    cat(" with correction:", object$Correction)
  }
  if (!is.null(object$Tree)) {
    cat("\nPhylogenetic or functional entropy was calculated according to the tree", object$Tree, "\n", fill=TRUE)
    cat("Entropy is", ifelse(object$Normalized, "normalized", "not normalized"), fill=TRUE)
  }
  cat("\nEntropy equals:", object$Total)
  return(invisible(NULL))
}
