CommunityProfile <-
function(FUN, NorP, q.seq = seq(0, 2, 0.1), 
         NumberOfSimulations = 0, Alpha = 0.05, BootstrapMethod = "Chao2015", ..., CheckArguments = TRUE) 
{
  if (CheckArguments) {
    CheckentropartArguments()
  }

  # Estimated profile
  Values <- vapply(q.seq, function(q) FUN(NorP, q, ..., CheckArguments = FALSE), 0)
  
  if (NumberOfSimulations > 0) {
    NsInt <- round(NorP)
    if (any(abs(NsInt-NorP) > sum(NorP)*.Machine$double.eps)) warning("Evaluation of the confidence interval of community profiles requires integer abundances in argument NorP. Abundances have been rounded.")
    # Create a MetaCommunity made of simulated communities
    MCSim <- rCommunity(NumberOfSimulations, NorP=NsInt, BootstrapMethod=BootstrapMethod, CheckArguments = FALSE)
    # May return NA if the bootstrap method is not recognized
    if (any(is.na(MCSim))) stop("Communities could not be simulated.")
    ProgressBar <- utils::txtProgressBar(min=0, max=NumberOfSimulations)
    Sims <- matrix(nrow=NumberOfSimulations, ncol=length(q.seq))
    # Loops are required for the progress bar, instead of:
    # Sims <- apply(MCSim$Nsi, 2, function(Nsi) CommunityProfile(FUN, Nsi, q.seq, ...)$y)
    for (i in 1:NumberOfSimulations) {
      # Parralelize. Do not allow more forks in PhyloApply()
      ProfileAsaList <- parallel::mclapply(q.seq, function(q) FUN(MCSim$Nsi[, i], q, ..., CheckArguments=FALSE), mc.allow.recursive=FALSE)
      Sims[i, ] <- simplify2array(ProfileAsaList)
      utils::setTxtProgressBar(ProgressBar, i)
    }
    # Recenter simulated values
    Means <- apply(Sims, 2, mean)
    Sims <- t(t(Sims)-Means+Values)
    
    # Quantiles of simulations for each q
    EstEnvelope <- apply(Sims, 2, stats::quantile, probs = c(Alpha/2, 1-Alpha/2))
    colnames(EstEnvelope) <- q.seq
    Profile <- list(x=q.seq,
                    y=Values,
                    low=EstEnvelope[1,],
                    high=EstEnvelope[2,])
  } else {
    Profile <- list(x=q.seq,
                    y=Values)
  }
  
  class(Profile) <- "CommunityProfile"
  return (Profile)
}


as.CommunityProfile <-
function (x, y, low = NULL, high = NULL) 
{
  if (!is.numeric(x))
    stop("x must be a numeric vector")
  if (!is.numeric(y))
    stop("y must be a numeric vector")
  if (length(x) != length(y))
    stop("x and y must have the same length")
  
  Profile <- list(x=x, y=y)
  if (!is.null(low)) {
    if (length(x) != length(low))
      stop("x and low must have the same length")
    Profile$low <- low
  }
  if (!is.null(high)) {
    if (length(x) != length(high))
      stop("x and high must have the same length")
    Profile$high <- high
  }
  class(Profile) <- "CommunityProfile"
  return(Profile)
}


is.CommunityProfile <-
function (x) 
{
  inherits(x, "CommunityProfile")
}


plot.CommunityProfile <- 
function(x, ..., main = NULL, 
         xlab = "Order of Diversity", ylab = "Diversity", ylim = NULL,
         LineWidth = 2, ShadeColor = "grey75", BorderColor = "red")
{  
  if (is.null(ylim)) {
    # Evaluate ylim if not set by an argument
    if (is.null(x$low)) {
      ymin <- min(x$y)
    } else {
      ymin <- min(x$low)
    }
    if (is.null(x$high)) {
      ymax <- max(x$y)
    } else {
      ymax <- max(x$high)
    }
  } else {
    ymin <- ylim[1]
    ymax <- ylim[2]
  }
  
  graphics::plot(x=x$x, y=x$y, type="n", main=main, xlab=xlab, ylab=ylab, ylim=c(ymin, ymax), ...)
  CEnvelope(x, LineWidth=LineWidth, ShadeColor=ShadeColor, BorderColor=BorderColor)
}

