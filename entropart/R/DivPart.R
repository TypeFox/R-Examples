DivPart <-
function(q = 1, MC, Biased = TRUE, Correction = "Best", Tree = NULL, Normalize = TRUE, Z = NULL, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Preprocess the tree
  ppTree <- Preprocess.Tree(Tree)
  if (Normalize) {
    Height <- 1
  } else {
    Height <- ppTree$Height
  }  

  # Alpha and beta entropy of communities
  if (Biased) {
    AlphaEntropy <- AlphaEntropy(MC, q, "None", ppTree, Normalize, Z, CheckArguments=FALSE)
    GammaEntropy <- GammaEntropy(MC, q, "None", ppTree, Normalize, Z, CheckArguments=FALSE)
    BetaEntropy  <- BetaEntropy (MC, q, "None", ppTree, Normalize, Z, CheckArguments=FALSE)
  } else {
    AlphaEntropy <- AlphaEntropy(MC, q, Correction, ppTree, Normalize, Z, CheckArguments=FALSE)
    GammaEntropy <- GammaEntropy(MC, q, Correction, ppTree, Normalize, Z, CheckArguments=FALSE)
    # beta is calculated as gamma-alpha to ensure continuity. Community beta entropy is not calculated.
    BetaEntropy  <- list(Communities = NA, Total = GammaEntropy - AlphaEntropy$Total)      
  }
  # Total Diversities
  AlphaDiversity <- expq(AlphaEntropy$Total / Height, q) * Height
  GammaDiversity <- expq(GammaEntropy / Height, q) * Height
  BetaDiversity  <- GammaDiversity / AlphaDiversity
    # equals: expq(BetaEntropy$Total / Height / (1 - (q-1)*AlphaEntropy$Total/Height), q) 
  
  DivPart <- (list(
    MetaCommunity = deparse(substitute(MC)),
    Order = q, 
    Biased = Biased, 
    Correction = Correction,
    Normalized = Normalize,
    TotalAlphaDiversity = AlphaDiversity, 
    TotalBetaDiversity = BetaDiversity, 
    GammaDiversity = GammaDiversity, 
    CommunityAlphaDiversities = expq(AlphaEntropy$Communities / Height, q) * Height, 
    TotalAlphaEntropy = AlphaEntropy$Total, 
    TotalBetaEntropy = BetaEntropy$Total, 
    GammaEntropy = GammaEntropy, 
    CommunityAlphaEntropies = AlphaEntropy$Communities, 
    CommunityBetaEntropies = BetaEntropy$Communities
    ))
  if(!is.null(Tree))
    DivPart$Tree <- deparse(substitute(Tree)) 
  if(is.null(Z)) {
    DivPart$Method <- "HCDT"
  } else {
    DivPart$Method <- "Similarity-based"
    DivPart$Z <- deparse(substitute(Z))  
  }
  class(DivPart) <- "DivPart"
  
  return (DivPart)
}


is.DivPart <-
function (x) 
{
  inherits(x, "DivPart")
}


plot.DivPart <- 
function (x, ...) 
{
  graphics::plot(c(0, x$GammaDiversity), c(0, length(x$CommunityAlphaDiversities)), type = "n", xlab = expression(paste(alpha, " and ", gamma, " diversity")), ylab = expression(paste(beta, " diversity")), ...)
  graphics::rect(0, 0, x$GammaDiversity, 1, lty=2)
  graphics::rect(0, 0, x$TotalAlphaDiversity, x$TotalBetaDiversity, lty=2)
}


summary.DivPart <-
function(object, ...) 
{    
  cat(object$Method, "diversity partitioning of order", object$Order, "of metaCommunity", object$MetaCommunity, fill=TRUE)
  if (!object$Biased)  
    cat(" with correction:", object$Correction)
  cat("\n")
  
  if (!is.null(object$Tree)) {
    cat("Phylogenetic or functional diversity was calculated\naccording to the tree", object$Tree, "\n", fill=TRUE)
    cat("Diversity is", ifelse(object$Normalized, "normalized", "not normalized"), "\n", fill=TRUE)
  }
  if (!is.null(object$Z)) {
    cat("Phylogenetic or functional entropy was calculated\naccording to the similarity matrix", object$Z, "\n", fill=TRUE)
  }
  
  cat("Alpha diversity of communities:", "\n")
  print(object$CommunityAlphaDiversities)
  cat("Total alpha diversity of the communities:", "\n")
  print(object$TotalAlphaDiversity)
  cat("Beta diversity of the communities:", "\n")
  print(object$TotalBetaDiversity)
  cat("Gamma diversity of the metacommunity:", "\n")
  print(object$GammaDiversity)
  
  return(invisible(NULL))
}