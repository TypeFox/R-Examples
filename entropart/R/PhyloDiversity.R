PhyloDiversity <-
function(NorP, q = 1, Tree, Normalize = TRUE, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  UseMethod("PhyloDiversity")
}


PhyloDiversity.ProbaVector <-
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
  
  ppTree <- Preprocess.Tree(Tree)
  if (Normalize) {
    Height <- 1
  } else {
    Height <- ppTree$Height
  }  
  # Calculate entropy
  Diversity <- PhyloEntropy(NorP, q, ppTree, Normalize=TRUE, CheckArguments=FALSE)
  # Transform it into diversity
  Diversity$Cuts <- expq(Diversity$Cuts, q)
  Diversity$Total <- expq(Diversity$Total, q) * Height
  # Complete it
  Diversity$Function <- "PhyloDiversity" 
  Diversity$Distribution <- deparse(substitute(NorP))
  Diversity$Tree <- deparse(substitute(Tree))
  Diversity$Type <- "alpha or gamma"
  Diversity$Order <- q
  Diversity$Normalized <- Normalize
  
  class(Diversity) <- c("PhyloDiversity", "PhyloValue")
  
  return(Diversity)  
}


PhyloDiversity.AbdVector <-
function(NorP, q = 1, Tree, Normalize = TRUE, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return(bcPhyloDiversity(Ns=NorP, q=q, Tree=Tree, Normalize=Normalize, Correction=Correction, CheckArguments=CheckArguments))
}


PhyloDiversity.integer <-
function(NorP, q = 1, Tree, Normalize = TRUE, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return(bcPhyloDiversity(Ns=NorP, q=q, Tree=Tree, Normalize=Normalize, Correction=Correction, CheckArguments=CheckArguments))
}


PhyloDiversity.numeric <-
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
    return(PhyloDiversity.ProbaVector(NorP, q=q, Tree=Tree, Normalize=Normalize, CheckArguments=CheckArguments))
  } else {
    # Abundances
    return(PhyloDiversity.AbdVector(NorP, q=q, Tree=Tree, Normalize=Normalize, Correction=Correction, CheckArguments=CheckArguments))
  }
}


bcPhyloDiversity <-
function(Ns, q = 1, Tree, Normalize = TRUE, Correction = "Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  ppTree <- Preprocess.Tree(Tree)
  if (Normalize) {
    Height <- 1
  } else {
    Height <- ppTree$Height
  }  
  
  # Calculate normalized entropy (Height will be considered later)
  Diversity <- bcPhyloEntropy(Ns, q, ppTree, Normalize=TRUE, Correction, CheckArguments=FALSE)
  # Transform it into diversity
  Diversity$Cuts <- expq(Diversity$Cuts, q)
  Diversity$Total <- expq(Diversity$Total, q) * Height
  # Complete it
  Diversity$Function <- "bcPhyloDiversity" 
  Diversity$Distribution <- deparse(substitute(Ns))
  Diversity$Tree <- deparse(substitute(Tree))
  Diversity$Type <- "alpha or gamma"
  Diversity$Order <- q
  Diversity$Correction <- Correction
  Diversity$Normalized <- Normalize
  
  class(Diversity) <- c("PhyloDiversity", "PhyloValue")
  
  return(Diversity)  
}


is.PhyloDiversity <-
function (x) 
{
  inherits(x, "PhyloDiversity")
}


summary.PhyloDiversity <-
function(object, ...) 
{
  cat(object$Type, "phylogenetic or functional diversity of order", object$Order, "of distribution", object$Distribution, fill=TRUE)
  if (!is.null(object$Correction)) {
    cat(" with correction:", object$Correction)
  }
  if (!is.null(object$Tree)) {
    cat("\nPhylogenetic or functional diversity was calculated according to the tree", object$Tree, "\n", fill=TRUE)
    cat("Diversity is", ifelse(object$Normalized, "normalized", "not normalized"), fill=TRUE)
  }
  cat("\nDiversity equals:", object$Total)
  return(invisible(NULL))
}