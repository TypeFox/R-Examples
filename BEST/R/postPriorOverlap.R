
postPriorOverlap <-
function( paramSampleVec, prior, ..., yaxt="n", ylab="",
           xlab="Parameter", main="", cex.lab=1.5, cex=1.4,
           xlim=range(paramSampleVec), breaks=NULL) {

  # Does a posterior histogram for a single parameter, adds the prior,
  #   displays and calculates the overlap.
  # Returns the overlap.

  oldpar <- par(xpd=NA) ; on.exit(par(oldpar))

  # get breaks: a sensible number over the hdi; cover the full range (and no more);
  #   equal spacing.
  if (is.null(breaks)) {
    nbreaks <- ceiling(diff(range(paramSampleVec)) / as.numeric(diff(hdi(paramSampleVec))/18))
    breaks <- seq(from=min(paramSampleVec), to=max(paramSampleVec), length.out=nbreaks)
  }
  # plot posterior histogram.
  histinfo <- hist(paramSampleVec, xlab=xlab, yaxt=yaxt, ylab=ylab,
                   freq=FALSE, border='white', col='skyblue',
                   xlim=xlim, main=main, cex=cex, cex.lab=cex.lab,
                   breaks=breaks)

  if (is.numeric(prior))  {
    # plot the prior if it's numeric
    priorInfo <- hist(prior, breaks=c(-Inf, breaks, Inf), add=TRUE,
      freq=FALSE, col='yellow', border='white')$density[2:length(breaks)]
  } else if (is.function(prior)) {
    if(class(try(prior(0.5, ...), TRUE)) == "try-error")
      stop(paste("Incorrect arguments for the density function", substitute(prior)))
    priorInfo <- prior(histinfo$mids, ...)
  }
  # get (and plot) the overlap
  minHt <- pmin(priorInfo, histinfo$density)
  rect(breaks[-length(breaks)], rep(0, length(breaks)-1), breaks[-1], minHt, col='green',
    border='white')
  overlap <- sum(minHt * diff(histinfo$breaks))
  # Add curve if prior is a function
  if (is.function(prior))
    lines(histinfo$mids, priorInfo, lwd=2, col='brown')
  # Add text
  text(mean(breaks), 0, paste0("overlap = ", round(overlap*100), "%"), pos=3, cex=cex)
    
  return(overlap)
}
