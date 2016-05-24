DivProfile <-
function(q.seq = seq(0, 2, .1), MC, Biased = TRUE, Correction = "Best", Tree = NULL, Normalize = TRUE, Z = NULL, NumberOfSimulations = 0, Alpha = 0.05, CheckArguments = TRUE) 
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

  # Calculate diversity profile. Parallelize.
  Diversity.seq <- simplify2array(mclapply(q.seq, DivPart, MC=MC, Biased=Biased, Correction=Correction, Tree=ppTree, Normalize=Normalize, Z=Z, CheckArguments =FALSE))

  # Rearrange complex structures
  Dalpha <- unlist(Diversity.seq["CommunityAlphaDiversities", ])
  arrDalpha <- simplify2array(tapply(Dalpha, names(Dalpha), c))
  row.names(arrDalpha) <- q.seq
  Ealpha <- unlist(Diversity.seq["CommunityAlphaEntropies", ])
  arrEalpha <- simplify2array(tapply(Ealpha, names(Ealpha), c))
  row.names(arrEalpha) <- q.seq
  # Prepare a list of results
  DivProfile <- list(MetaCommunity = unlist(Diversity.seq["MetaCommunity", 1], use.names=FALSE),
                     Order = unlist(Diversity.seq["Order", ], use.names=FALSE), 
                     Biased = unlist(Diversity.seq["Biased", 1], use.names=FALSE), 
                     Correction = unlist(Diversity.seq["Correction", 1], use.names=FALSE),
                     Normalized = unlist(Diversity.seq["Normalized", 1], use.names=FALSE),
                     CommunityAlphaDiversities = arrDalpha, 
                     CommunityAlphaEntropies = arrEalpha, 
                     TotalAlphaDiversity = unlist(Diversity.seq["TotalAlphaDiversity", ]), 
                     TotalBetaDiversity = unlist(Diversity.seq["TotalBetaDiversity", ]), 
                     GammaDiversity = unlist(Diversity.seq["GammaDiversity", ]), 
                     TotalAlphaEntropy =  unlist(Diversity.seq["TotalAlphaEntropy", ]), 
                     TotalBetaEntropy = unlist(Diversity.seq["TotalBetaEntropy", ]), 
                     GammaEntropy =  unlist(Diversity.seq["GammaEntropy", ])
                    )
  if(!is.null(Tree))
    DivProfile$Tree <- deparse(substitute(Tree)) 
  if(is.null(Z)) {
    DivProfile$Method <- "HCDT"
  } else {
    DivProfile$Method <- "Similarity-based"
    DivProfile$Z <- deparse(substitute(Z))  
  }
  class(DivProfile) <- "DivProfile"
  
  # Confidence interval
  if (NumberOfSimulations > 0) {
    
    Q <- length(q.seq)
    ProgressBar <- utils::txtProgressBar(min=-3, max=Q, initial=-3)
    ## Obtain simulated Metacommunities
    RedrawSpecies<- function(SpeciesAbundances){
      # Very simplified (for speed) version of rCommunity with BootstrapMethod="Marcon"
      MetaCommunity(stats::rmultinom(NumberOfSimulations, sum(SpeciesAbundances), SpeciesAbundances))
    }
    utils::setTxtProgressBar(ProgressBar, -2)
    # Resample each community according to species abundances
    ResampledCs <- apply(MC$Nsi, 2, function(Ns) RedrawSpecies(SpeciesAbundances=Ns))
    utils::setTxtProgressBar(ProgressBar, -1)
    # Each MC of this list is a simulation set of each original community
    # Build simulated MCs by picking simulated communities
    SimMC <- lapply(1:NumberOfSimulations, function(i) MetaCommunity(sapply(ResampledCs, function(mc) mc$Nsi[, i]), Weights=MC$Wi))
    utils::setTxtProgressBar(ProgressBar, -0)
    
    ## Calculate alpha, beta and gamma of each simulated MC, at each q
    # Prepare a matrix to store envelopes.
    Envelopes <- matrix(nrow=12, ncol=Q)
    rownames(Envelopes) <- c("TotalAlphaEntropyLow", "TotalAlphaEntropyHigh", 
                             "TotalBetaEntropyLow", "TotalBetaEntropyHigh",
                             "GammaEntropyLow", "GammaEntropyHigh",
                             "TotalAlphaDiversityLow", "TotalAlphaDiversityHigh", 
                             "TotalBetaDiversityLow", "TotalBetaDiversityHigh",
                             "GammaDiversityLow", "GammaDiversityHigh")
    for (qi in 1:Q) {
      # Calculate diversites. Parallelize.
      Diversity.qi <- simplify2array(mclapply(SimMC, function(mc) DivPart(q=q.seq[qi], MC=mc, Biased=Biased, Correction=Correction, Tree=ppTree, Normalize=Normalize, Z=Z, CheckArguments=FALSE)))
      # Put alpha and gamma simulated entropies into vectors
      TotalAlphaEntropy <-  unlist(Diversity.qi["TotalAlphaEntropy", ])
      TotalBetaEntropy <-  unlist(Diversity.qi["TotalBetaEntropy", ])
      GammaEntropy <-  unlist(Diversity.qi["GammaEntropy", ])
      # Recenter entropy
      TotalAlphaEntropy <- TotalAlphaEntropy + DivProfile$TotalAlphaEntropy[qi] - mean(TotalAlphaEntropy, na.rm=TRUE)
      TotalBetaEntropy <- TotalBetaEntropy + DivProfile$TotalBetaEntropy[qi] - mean(TotalBetaEntropy, na.rm=TRUE)
      GammaEntropy <- GammaEntropy + DivProfile$GammaEntropy[qi] - mean(GammaEntropy, na.rm=TRUE)
      TotalAlphaDiversity <- expq(TotalAlphaEntropy, q.seq[qi])
      GammaDiversity <- expq(GammaEntropy, q.seq[qi])
      TotalBetaDiversity <- GammaDiversity / TotalAlphaDiversity
      # Quantiles
      Envelopes["TotalAlphaEntropyLow", qi] <- stats::quantile(TotalAlphaEntropy, probs=Alpha, na.rm=TRUE)
      Envelopes["TotalAlphaEntropyHigh", qi] <- stats::quantile(TotalAlphaEntropy, probs=1-Alpha, na.rm=TRUE)
      Envelopes["TotalBetaEntropyLow", qi] <- stats::quantile(TotalBetaEntropy, probs=Alpha, na.rm=TRUE)
      Envelopes["TotalBetaEntropyHigh", qi] <- stats::quantile(TotalBetaEntropy, probs=1-Alpha, na.rm=TRUE)
      Envelopes["GammaEntropyLow", qi] <- stats::quantile(GammaEntropy, probs=Alpha, na.rm=TRUE)
      Envelopes["GammaEntropyHigh", qi] <- stats::quantile(GammaEntropy, probs=1-Alpha, na.rm=TRUE)
      Envelopes["TotalAlphaDiversityLow", qi] <- stats::quantile(TotalAlphaDiversity, probs=Alpha, na.rm=TRUE)
      Envelopes["TotalAlphaDiversityHigh", qi] <- stats::quantile(TotalAlphaDiversity, probs=1-Alpha, na.rm=TRUE)
      Envelopes["TotalBetaDiversityLow", qi] <- stats::quantile(TotalBetaDiversity, probs=Alpha, na.rm=TRUE)
      Envelopes["TotalBetaDiversityHigh", qi] <- stats::quantile(TotalBetaDiversity, probs=1-Alpha, na.rm=TRUE)
      Envelopes["GammaDiversityLow", qi] <- stats::quantile(GammaDiversity, probs=Alpha, na.rm=TRUE)
      Envelopes["GammaDiversityHigh", qi] <- stats::quantile(GammaDiversity, probs=1-Alpha, na.rm=TRUE)
      # Progressbar
      utils::setTxtProgressBar(ProgressBar, qi)
    }
    
    # Integrate the envelopes into the object
    DivProfile$TotalAlphaEntropyLow <- Envelopes["TotalAlphaEntropyLow", ]
    DivProfile$TotalAlphaEntropyHigh <- Envelopes["TotalAlphaEntropyHigh", ]
    DivProfile$TotalBetaEntropyLow <- Envelopes["TotalBetaEntropyLow", ]
    DivProfile$TotalBetaEntropyHigh <- Envelopes["TotalBetaEntropyHigh", ]
    DivProfile$GammaEntropyLow <- Envelopes["GammaEntropyLow", ]
    DivProfile$GammaEntropyHigh <- Envelopes["GammaEntropyHigh", ]
    DivProfile$TotalAlphaDiversityLow <- Envelopes["TotalAlphaDiversityLow", ]
    DivProfile$TotalAlphaDiversityHigh <- Envelopes["TotalAlphaDiversityHigh", ]
    DivProfile$TotalBetaDiversityLow <- Envelopes["TotalBetaDiversityLow", ]
    DivProfile$TotalBetaDiversityHigh <- Envelopes["TotalBetaDiversityHigh", ]
    DivProfile$GammaDiversityLow <- Envelopes["GammaDiversityLow", ]
    DivProfile$GammaDiversityHigh <- Envelopes["GammaDiversityHigh", ]
  }
  
  return (DivProfile)
}


is.DivProfile <-
function (x) 
{
  inherits(x, "DivProfile")
}


plot.DivProfile <- 
function (x, ..., main = NULL, xlab = "Order of Diversity", ylab = NULL, Which = "All", 
          LineWidth = 2, ShadeColor = "grey75", BorderColor = "red")
{
  # Save graphical parameters
  if (Which == "All") {
    op <- graphics::par(no.readonly = TRUE)
    graphics::par(mfrow=c(2, 2))    
  }
  if (Which == "All" | (Which == "Alpha" & is.null(main))) main <- "Total Alpha Diversity"
  if (Which == "All" | (Which == "Alpha" & is.null(ylab))) ylab <- expression(paste(alpha, " diversity"))
  if (Which == "All" | Which == "Alpha") {
    graphics::plot(y=x$TotalAlphaDiversity, x=x$Order, type="n", main=main, xlab=xlab, ylab=ylab, ylim=c(min(x$TotalAlphaDiversity, x$TotalAlphaDiversityLow, na.rm=TRUE), max(x$TotalAlphaDiversity, x$TotalAlphaDiversityHigh, na.rm=TRUE)), ...)
    if (!(is.null(x$TotalAlphaDiversityHigh) | is.null(x$TotalAlphaDiversityLow))) {
      # Shaded polygon (adapted from Didzis Elferts, 
      # http://stackoverflow.com/questions/14069629/plotting-confidence-intervals)
      graphics::polygon(c(x$Order, rev(x$Order)), c(pmax(x$TotalAlphaDiversityLow, graphics::par('usr')[3]), pmin(rev(x$TotalAlphaDiversityHigh), par('usr')[4])), col = ShadeColor, border = FALSE)
      # Add red lines on borders of polygon
      graphics::lines(x$Order, x$TotalAlphaDiversityHigh, col=BorderColor, lty=2)
      graphics::lines(x$Order, x$TotalAlphaDiversityLow, col=BorderColor, lty=2)
    }
    graphics::lines(y=x$TotalAlphaDiversity, x=x$Order, main=main, xlab=xlab, ylab=ylab, lwd=LineWidth, ...)
  }
  if (Which == "All" | (Which == "Communities" & is.null(main))) main <- "Alpha Diversity of Communities"
  if (Which == "All" | (Which == "Communities" & is.null(ylab))) ylab <- expression(paste(alpha, " diversity"))
  if (Which == "All" | Which == "Communities") {
    Palette <- grDevices::palette(grDevices::rainbow(ncol(x$CommunityAlphaDiversities)))
    graphics::plot(x$CommunityAlphaDiversities[, 1] ~ x$Order, type="n", xlim=c(min(x$Order), max(x$Order)), ylim=c(min(x$CommunityAlphaDiversities), max(x$CommunityAlphaDiversities)), main=main, xlab=xlab, ylab=ylab, ...)
    for (Community in (1:ncol(x$CommunityAlphaDiversities))) {
      graphics::lines(x=x$Order, y=x$CommunityAlphaDiversities[, Community], lty=Community, col=Palette[Community])
    }  
    if (Which == "Communities") {
      graphics::legend("topright", colnames(x$CommunityAlphaDiversities), lty=1:ncol(x$CommunityAlphaDiversities), col=Palette, inset=0.01)
    }
  }
  if (Which == "All" | (Which == "Beta" & is.null(main))) main <- "Beta Diversity"
  if (Which == "All" | (Which == "Beta" & is.null(ylab))) ylab <- expression(paste(beta, " diversity"))
  if (Which == "All" | Which == "Beta") {
    graphics::plot(y=x$TotalBetaDiversity, x=x$Order, type="n", main=main, xlab=xlab, ylab=ylab, ylim=c(min(x$TotalBetaDiversity, x$TotalBetaDiversityLow, na.rm=TRUE), max(x$TotalBetaDiversity, x$TotalBetaDiversityHigh, na.rm=TRUE)), ...)
    if (!(is.null(x$TotalBetaDiversityHigh) | is.null(x$TotalBetaDiversityLow))) {
      # Shaded polygon (adapted from Didzis Elferts, 
      # http://stackoverflow.com/questions/14069629/plotting-confidence-intervals)
      graphics::polygon(c(x$Order, rev(x$Order)), c(pmax(x$TotalBetaDiversityLow, graphics::par('usr')[3]), pmin(rev(x$TotalBetaDiversityHigh), par('usr')[4])), col = ShadeColor, border = FALSE)
      # Add red lines on borders of polygon
      graphics::lines(x$Order, x$TotalBetaDiversityHigh, col=BorderColor, lty=2)
      graphics::lines(x$Order, x$TotalBetaDiversityLow, col=BorderColor, lty=2)
    }
    graphics::lines(y=x$TotalBetaDiversity, x=x$Order, main=main, xlab=xlab, ylab=ylab, lwd=LineWidth, ...)
  }
  if (Which == "All" | (Which == "Gamma" & is.null(main))) main <- "Gamma Diversity"
  if (Which == "All" | (Which == "Gamma" & is.null(ylab))) ylab <- expression(paste(gamma, " diversity"))
  if (Which == "All" | Which == "Gamma") {
    graphics::plot(y=x$GammaDiversity, x=x$Order, type="n", main=main, xlab=xlab, ylab=ylab, ylim=c(min(x$GammaDiversity, x$GammaDiversityLow, na.rm=TRUE), max(x$GammaDiversity, x$GammaDiversityHigh, na.rm=TRUE)), ...)
    if (!(is.null(x$GammaDiversityHigh) | is.null(x$GammaDiversityLow))) {
      # Shaded polygon (adapted from Didzis Elferts, 
      # http://stackoverflow.com/questions/14069629/plotting-confidence-intervals)
      graphics::polygon(c(x$Order, rev(x$Order)), c(pmax(x$GammaDiversityLow, graphics::par('usr')[3]), pmin(rev(x$GammaDiversityHigh), par('usr')[4])), col = ShadeColor, border = FALSE)
      # Add red lines on borders of polygon
      graphics::lines(x$Order, x$GammaDiversityHigh, col=BorderColor, lty=2)
      graphics::lines(x$Order, x$GammaDiversityLow, col=BorderColor, lty=2)
    }
    graphics::lines(y=x$GammaDiversity, x=x$Order, main=main, xlab=xlab, ylab=ylab, lwd=LineWidth, ...)
  }
  # Restore parameters
  if (Which == "All") {
    graphics::par(op)
  }
}


summary.DivProfile <-
function(object, ...) 
{
  cat("Diversity profile of MetaCommunity", object$MetaCommunity, fill=TRUE)
  if (!object$Biased)  
    cat(" with correction:", object$Correction)
  cat("\n")
  
  if (!is.null(object$Tree)) {
    cat("Phylogenetic or functional diversity was calculated according to the tree", object$Tree, "\n", fill=TRUE)
    cat("Diversity is", ifelse(object$Normalized, "normalized", "not normalized"), "\n", fill=TRUE)
  }
  
  cat("Diversity against its order:\n")
  Values <- cbind(object$Order, object$TotalAlphaDiversity, object$TotalBetaDiversity, object$GammaDiversity)
  colnames(Values) <- c("Order", "Alpha Diversity", "Beta Diversity", "Gamma Diversity")
  print(Values)
  
  return(invisible(NULL))
}