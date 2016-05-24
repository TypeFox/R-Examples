plot.copas <- function(x,
                       which=1:4,
                       caption=c("Funnel plot", "Contour plot",
                         "Treatment effect plot",
                         "P-value for residual selection bias"),
                       xlim.pp=NULL,
                       level=0.95,
                       orthogonal.line=TRUE,
                       lines=FALSE,
                       sign.rsb=x$sign.rsb,
                       warn=-1,
                       ...){
  
  meta:::chkclass(x, "copas")
  
  if (!is.numeric(which) || any(which < 1) || any(which > 4)) 
    stop("'Argument which' must be in 1:4")
  
  if (is.null(sign.rsb))
    sign.rsb <- 0.1
  else
    meta:::chklevel(sign.rsb)
  
  
  oldwarn <- options()$warn
  on.exit(options(warn=oldwarn))
  options(warn=warn)
  
  
  TE <- x$TE
  seTE <- x$seTE
  TE.random <- x$TE.random
  seTE.random <- x$seTE.random
  gamma0 <- x$gamma0
  gamma1 <- x$gamma1
  TE.contour <- x$TE.contour
  levels <- x$regr$levels
  slope <- x$slope
  x.slope <- x$x.slope
  y.slope <- x$y.slope
  TE.slope <- x$TE.slope
  seTE.slope <- x$seTE.slope
  publprob <- x$publprob
  pval.rsb <- x$pval.rsb
  sm <- x$sm
  ##
  ord <- order(publprob)
  
  
  show <- rep(FALSE, 4)
  show[which] <- TRUE
  ##
  if (sum(show) == 1) oldpar <- par(pty="s")
  if (sum(show) == 2) oldpar <- par(mfrow=c(1,2), pty="s")
  if (sum(show) >  2) oldpar <- par(mfrow=c(2,2),
                                       mar=c(4.0,4.1,1.5,0.5),
                                       pty="s")
  ##
  on.exit(par(oldpar), add=TRUE)
  
  
  if (meta:::is.relative.effect(sm))
    sm <- paste("log ", sm, sep="")
  
  
  ## Plot 1: funnel plot
  ##
  if (show[1]){
    funnel(x=TE, y=seTE, level=0.95, comb.fixed=TRUE, sm=sm, xlab="", ylab="")
    abline(v=TE.random, lty=1, col="darkgray")
    ##
    mtext(sm, side=1, line=2)
    mtext("Standard error", side=2, line=2)
    box()
    title(main=caption[1])
  }
  
  
  ## Plot 2: contour plot
  ##
  ## NB this is a square plot, with axes transformed.
  ##
  gamma0.min <- min(gamma0)
  gamma0.max <- max(gamma0)
  gamma1.min <- min(gamma1)
  gamma1.max <- max(gamma1)
  ##
  gamma0.rescale <- (gamma0-gamma0.min) / (gamma0.max-gamma0.min)
  gamma1.rescale <- (gamma1-gamma1.min) / (gamma1.max-gamma1.min)
  ##
  if (show[2]){
    contour(gamma0.rescale, gamma1.rescale, TE.contour,
            xlim=c(0, 1), ylim=c(0, 1),
            xlab="", ylab="",
            labcex=0.8, lty=2, col="black", axes=FALSE, cex=2,
            levels=levels)
    ##
    axis(side=1, at=axTicks(1),
         labels=round(axTicks(1)*(gamma0.max-gamma0.min)+gamma0.min, 2))
    axis(side=2, at=axTicks(2),
         labels=round(axTicks(2)*(gamma1.max-gamma1.min)+gamma1.min, 2))
    ##
    mtext("Values of gamma0", side=1, line=2)
    mtext("Values of gamma1", side=2, line=2)
    ##
    box()
    ##
    ## text(1,1, round(TE.random, 2)) not particularly useful
    ##
    ## Add estimated orthogonal line
    ##
    if (orthogonal.line){
      lines(gamma0.rescale, 1 + slope*(gamma0.rescale-1), lty=1)
      points(x.slope, y.slope)
    }
    ##
    if (lines){
      if (is.null(x$min.r.squared))
        min.r.squared <- -100
      else
        min.r.squared <- x$min.r.squared
      ##
      for (i in seq(along=x$regr$intercepts)[x$regr$adj.r.squareds>0])
        abline(x$regr$intercepts[i], x$regr$slopes[i], col="green")
      for (i in seq(along=x$regr$intercepts)[x$regr$adj.r.squareds<=0])
        abline(x$regr$intercepts[i], x$regr$slopes[i], col="red", lty=3)
    }
    ##
    title(main=caption[2])
  }
  
  
  ## Plot 3:
  ## mean (plus level%-CI) against prob of publishing trial
  ## with largest SD down orthogonal line
  ##
  xvalue <- publprob[ord]
  yvalue <- TE.slope[ord]
  ##
  ci.y <- ci(yvalue, seTE.slope[ord], level=level)
  ci.y$lower[is.infinite(ci.y$lower)] <- NA
  ci.y$upper[is.infinite(ci.y$upper)] <- NA
  ##
  if (show[3]){
    if (is.null(xlim.pp))
      xlim <- range(xvalue, na.rm=TRUE)[c(2,1)]
    else
      xlim <- xlim.pp
    ##
    xlim[is.infinite(xlim)] <- c(1,0)[is.infinite(xlim)]
    ##
    if (diff(xlim)==0)
      xlim <- xlim + c(0.0001, -0.0001)
    ##
    ##
    plot(xvalue, yvalue,
         type="l",
         xlim=xlim,
         ylim=c(
           min(c(0, ci.y$lower), na.rm=TRUE),
           max(c(ci.y$upper, 0), na.rm=TRUE)
           ),
         xlab="", ylab="", axes=FALSE)
    ##
    axis(side=1, axTicks(1), labels=round(axTicks(1), 2))
    axis(side=2, axTicks(2), labels=round(axTicks(2), 2))
    ##
    mtext("Probability of publishing the trial with largest sd",
          side=1, line=2)
    mtext(sm, side=2, line=2)
    ##
    points(xvalue, yvalue)
    ##
    abline(h=0)
    ##
    ci.random <- ci(TE.random, seTE.random, level=level)
    abline(h=TE.random, lty=1, col="darkgray")
    abline(h=ci.random$lower, lty=1, col="gray")
    abline(h=ci.random$upper, lty=1, col="gray")
    ##
    lines(xvalue, ci.y$lower, lty=1)
    lines(xvalue, ci.y$upper, lty=1)
    box()
    title(main=caption[3])
  }
  
  
  ## Plot 4:
  ## goodness of fit (p-value) against prob of publishing
  ## trial with largest SD
  ##
  yvalue <- pval.rsb[ord]
  sel.y <- !is.na(yvalue)
  xvalue <- xvalue[sel.y]
  yvalue <- yvalue[sel.y]
  ##
  if (show[4]){
    if (is.null(xlim.pp))
      xlim <- range(xvalue, na.rm=TRUE)[c(2,1)]
    else
      xlim <- xlim.pp
    ##
    xlim[is.infinite(xlim)] <- c(1,0)[is.infinite(xlim)]
    ##
    if (diff(xlim)==0)
      xlim <- xlim + c(0.0001, -0.0001)
    ##
    ##
    if (length(xvalue)==0){
      plot(xvalue, yvalue,
           xlim=xlim, ylim=c(0,1),
           xlab="", ylab="", axes=FALSE)
    }
    else{
      plot(loess(yvalue~xvalue),
           type="l",
           xlim=xlim, ylim=c(0,1),
           xlab="", ylab="", axes=FALSE)
    }
    ##
    points(xvalue, yvalue)
    ##
    axis(side=1, axTicks(1), labels=round(axTicks(1), 2))
    axis(side=2, axTicks(2), labels=round(axTicks(2), 2))
    ##
    mtext("Probability of publishing the trial with largest sd",
          side=1, line=2)
    mtext("P-value for residual selection bias", side=2, line=2)
    ##
    abline(h=sign.rsb, lty=2)
    box()
    ##
    title(main=caption[4])
  }
  
  
  invisible(NULL)
}
