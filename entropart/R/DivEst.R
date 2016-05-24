DivEst <-
function(q = 0, MC, Biased = TRUE, Correction = "Best", Tree = NULL, Normalize = TRUE, Z = NULL, Simulations = 100, CheckArguments = TRUE) 
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

  # Estimation from data
  RealEst <- DivPart(q, MC, Biased, Correction, ppTree, Normalize, Z, CheckArguments=FALSE)

  # RedrawSpecies resamples a community according to species abundances.
  RedrawSpecies<- function(SpeciesAbundances){
    # Very simplified (for speed) version of rCommunity with BootstrapMethod="Marcon"
    stats::rmultinom(1, sum(SpeciesAbundances), SpeciesAbundances)
  }

  # SimulateEntropy resamples all communities and calculates entropy
  SimulateEntropy <- function(Progression) {
    SimNsi <- apply(MC$Nsi, 2, RedrawSpecies)
    rownames(SimNsi) <- rownames(MC$Nsi) 
    SimMC <- Preprocess.MC(SimNsi, MC$Wi)
    NewSim <- DivPart(q, SimMC, Biased, Correction, Tree, Normalize, Z, CheckArguments=FALSE)
    # update progress bar
    utils::setTxtProgressBar(ProgressBar, Progression)
    c(NewSim$TotalAlphaEntropy, NewSim$TotalBetaEntropy, NewSim$GammaEntropy)
  }
  
  # Simulate entropy
  ProgressBar <- utils::txtProgressBar(min=0, max=Simulations)
  RawSimulatedEntropy <- sapply(1:Simulations, SimulateEntropy)
  # Recenter entropy
  SimulatedEntropy <- RawSimulatedEntropy+with(RealEst, c(TotalAlphaEntropy, TotalBetaEntropy, GammaEntropy))-apply(RawSimulatedEntropy, 1, mean)
  # Transform entropy to diversity
  if (q == 1) { 
    SimulatedDiversity <- exp(SimulatedEntropy)
  } else {
    SimulatedDiversity <- SimulatedEntropy
    SimulatedDiversity[1,] <- expq(SimulatedEntropy[1,] / Height, q) * Height
    SimulatedDiversity[3,] <- expq(SimulatedEntropy[3,] / Height, q) * Height
    SimulatedDiversity[2,] <- SimulatedDiversity[3,] / SimulatedDiversity[1,]
      # equals: expq(SimulatedEntropy[2,] / Height / (1 - (q-1)*SimulatedEntropy[1,]/Height), q) * Height
  }
  dimnames(SimulatedDiversity)[[1]] <- list("Alpha", "Beta", "Gamma")

  DivEst <- RealEst
  DivEst$SimulatedDiversity <- SimulatedDiversity
  class(DivEst) <- c("DivEst", class(RealEst))

  return (DivEst)
}


is.DivEst <-
function (x) 
{
  inherits(x, "DivEst")
}


plot.DivEst <- 
function (x, ..., main = NULL, Which = "All") 
{
  # Save graphical parameters
  op <- graphics::par(no.readonly = TRUE)
  graphics::par(mfrow=c(2, 2))
  
  if (Which == "All" | (Which == "Alpha" & is.null(main))) main <- "Alpha Diversity"
  if (Which == "All" | Which == "Alpha") {
    graphics::plot(as.SimTest(x$TotalAlphaDiversity, x$SimulatedDiversity["Alpha",]), main=main, ...)
  }
  if (Which == "All" | (Which == "Beta" & is.null(main))) main <- "Beta Diversity"
  if (Which == "All" | Which == "Beta") {
    graphics::plot(as.SimTest(x$TotalBetaDiversity, x$SimulatedDiversity["Beta",]), main=main, ...)
  }
  if (Which == "All" | (Which == "Gamma" & is.null(main))) main <- "Gamma Diversity"
  if (Which == "All" | Which == "Gamma") {
    graphics::plot(as.SimTest(x$GammaDiversity, x$SimulatedDiversity["Gamma",]), main=main, ...)
  }
  
  # Legend and restore parameters
  if (Which == "All") {
    graphics::par(mar=c(0, 0, 0, 0))
    graphics::plot(0:10, 0:10, type="n", xlab=NULL, frame.plot=FALSE, xaxt="n", yaxt="n", col.lab="white")
    leg <- c("Null Distribution", "True Estimate", "95% confidence interval") 
    graphics::legend(2, 8, leg, col = c(1, 2, 1), lty = 1:3, merge = TRUE, cex=1)
    graphics::par(op)
  }
}


summary.DivEst <-
function(object, ...) 
{
  cat("Diversity partitioning of order", object$Order, "of MetaCommunity", object$MetaCommunity, fill=TRUE)
  if (!object$Biased)  
    cat(" with correction:", object$Correction)
  cat("\n")
  
  if (!is.null(object$Tree)) {
    cat("Phylogenetic or functional diversity was calculated according to the tree", object$Tree, "\n", fill=TRUE)
    cat("Diversity is", ifelse(object$Normalized, "normalized", "not normalized"), "\n", fill=TRUE)
  }
  cat("Alpha diversity of communities:", "\n")
  print(object$CommunityAlphaDiversities)
  cat("Total alpha diversity of the communities:", "\n")
  print(object$TotalAlphaDiversity)
  cat("Beta diversity of the communities:", "\n")
  print(object$TotalBetaDiversity)
  cat("Gamma diversity of the metacommunity:", "\n")
  print(object$GammaDiversity)
  
  cat("Quantiles of simulations (alpha, beta and gamma diversity):\n")
  quant <- c(0, 0.01, 0.025, 0.05, 0.1, seq(0.25, 0.75, 0.25), 0.9, 0.95, 0.975, 0.99, 1)
  print(stats::quantile(object$SimulatedDiversity["Alpha", ], probs = quant))
  print(stats::quantile(object$SimulatedDiversity["Beta", ], probs = quant))
  print(stats::quantile(object$SimulatedDiversity["Gamma", ], probs = quant))
  
  return(invisible(NULL))
}