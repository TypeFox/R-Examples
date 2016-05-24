
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


############################################################################### 
# FUNCTION:              DESCRIPTION: 
#  stabilityAnalytics      Retroactive stability analytics 
# FUNCTION:              DESCRIPTION:        
#  turnsAnalytics          Retroactive turning point analytics 
#  drawdownsAnalytics      Retroactive maximum drawdown analytics
#  garchAnalytics          Retroactive Garch volatility analytics
#  riskmetricsAnalytics    Retroactive Riskmetrics analytics
#  bcpAnalytics            Retroactive Bayesian changepoints analytics
#  pcoutAnalytics          Retroactive Principal component outlier analytics
# FUNCTION:              DESCRIPTION:
#  addRainbow              Adds rainbow colored stability indicators
# FUNCTION:              DESCRIPTION: 
#  waveletSpectrum         Retroactive Morlet wavelet analytics
# FUNCTION:              DESCRIPTION:
#  parAnalytics            Graph frame settings for a desired analytics
############################################################################### 


stabilityAnalytics <- 
  function(index, method=c("turns", "drawdowns", "garch", 
                           "riskmetrics", "bcp", "pcout"), ...)
  {
    # A function implemented by Diethelm Wuertz and Tobias Setz
    
    # Description:
    # Retroactive stability analytics 
    
    # FUNCTION:
    
    # Run Analytics:
    method <- match.arg(method)
    FUN <- paste(method, "Analytics", sep = "")
    fun <- match.fun(FUN)
    
    # Return Value:
    fun(index, ...)
  }


###############################################################################


turnsAnalytics <- 
  function(index, spar=0.5, main=NULL, 
           trace=TRUE, doplot=TRUE, at=pretty(index), format="%m/%y")
  {
    # A function implemented by Diethelm Wuertz and Tobias Setz
    
    # Description:
    #   Retroactive Turning Points Analytics
    
    # Arguments:
    #   index - an index or price S4 'timeSeries' object
    #   spar - a numeric between 0 and, the degree of smoothness
    #   main - main plot title
    #   trace - a logical, should the results be traced?
    #   plot -  a logical, should the results be plotted
    #   at - generate pretty axis positions
    #   format -  a string describing the label format
    
    # Note:
    #   The lowess and supsmu smoothers are by far not as good as the
    #   spline smoother.
    
    # FUNCTION:
    
    # Settings:
    stopifnot(isUnivariate(index))
    if (is.null(main)) main <- "Retroactive Turnpoints Analytics"
    
    # Smooth Return Series:
    rets <- returns(index)
    indexSmu <- smoothSpline(log(index), spar=spar)
    
    # Turnpoints:
    warn <- getOption("warn")
    options(warn=-1)
    indexTurns <- turns(indexSmu[, 2])
    options(warn=warn)
    indexTurns <- indexTurns[indexTurns[, 2] !=0, ]
    
    # Add Verical Lines:
    turns.tps <- format(time(indexTurns))
    n.tps <- length(turns.tps)
    
    # Trace Results:
    if (trace){
      cat("Series:\n", colnames(index), "\n")
      cat("Turning Points:", n.tps, "\n")
      print(turns.tps) 
    }
    
    # Positions:
    positions <- sign(returns(exp(indexSmu[, 2]), trim=FALSE))
    positions[1, 1] <- 0
    ablines <- time(positions)[as.vector(positions) < 0]
    
    # Plot:
    if (doplot) {
      
      # Turnpoints:
      range <- range(log(index))
      ylim <- c(range[1], range[2] + diff(range(log(index)))/4)
      plot(indexSmu[, 1], ylim=ylim, main="", ylab="", xlab="", las=2, 
           at=at, format=format, col="black")
      title(main=main, ylab=paste("log Index", colnames(index)), xlab="")
      abline(v=ablines, lty=3, lwd=2, col="steelblue")
      lines(indexSmu[, 1], col="black") 
      lines(indexSmu[, 2], col="red") 
      if (turns.tps[2] > 0) points(indexTurns[, 1], pch=19, col="red")
      abline(v=as.POSIXct(at), col="darkgrey", lty=3)
      box(col="white")
      box(bty="l")
      
      # Add Returns:
      center <- range[2] + diff(range(log(index)))/4/2
      scale <- diff(range(log(index)))/4
      returnsScaled <- (rets-mean(rets))/max(abs(rets)) * scale/2 + center
      lines(returnsScaled, col="orange")
      abline(h=mean(returnsScaled), col="darkgrey", lty=3)
      box(col="white")
      box(bty="l")
      
      # Add Rmwtrics - Do not Remove!
      mtext("Rmetrics", adj=0, side=4, cex=0.7, col="darkgrey")
      
    }
    
    # Check Turnpoints:
    if ( turns.tps[1] == "1" & turns.tps[2] == "0")
      n.tps <- 0
    
    # Return Value:
    invisible(list(
      data=indexSmu, 
      turns=turns.tps, 
      positions=positions, 
      ablines=ablines, 
      n=n.tps, 
      smooth=spar))
  }


###############################################################################    


drawdownsAnalytics <- 
  function(index, spar=0.5, main=NULL, 
           trace=TRUE, doplot=TRUE, at=pretty(index), format="%m/%y")
  {
    # A function implemented by Diethelm Wuertz and Tobias Setz
    
    # Description:
    #   Retroactive Maximum Drawdown Analytics
    
    # Arguments:
    #   index - an index or price S4 'timeSeries' object
    #   spar - a numeric between 0 and, the degree of smoothness
    #   main - main plot title
    #   trace - a logical, should the results be traced?
    #   plot -  a logical, should the results be plotted
    #   at - generate pretty axis positions
    #   format -  a string describing the label format
    
    # FUNCTION:
    
    # Settings:
    stopifnot(isUnivariate(index))
    if (is.null(main)) main <- "Retroactive Drawdowns Analytics"
    
    # Series:
    rets<- returns(index)
    
    # Turning Points:
    tps <- turnsAnalytics(index=index, spar=spar, trace=trace, doplot=FALSE)
    ablines <- tps$ablines
    
    # Drawdowns:
    maxdd <- drawdowns(rets)
    
    # Plot:
    if(doplot) {
      
      # Plot:
      plot(maxdd, main=main, xlab="", ylab=paste("Drawdwons", colnames(index)),
           las=2, at=at, format=format)
      abline(v=ablines, lty=3, lwd=2, col="steelblue")
      lines(maxdd)
      box(col="white")
      box(bty="l")
      
      # Add Rmwtrics - Do not Remove!
      mtext("Rmetrics", adj=0, side=4, cex=0.7, col="darkgrey")
      
    }
    
    # Return Value:
    invisible(list(
      index = index, 
      series = maxdd))
  }


###############################################################################


garchAnalytics <-
  function (index, spar = 0.5, main=NULL, 
            trace=TRUE, doplot=TRUE, at=pretty(index), format="%m/%y")
  {
    # A function implemented by Diethelm Wuertz and Tobias Setz
    
    # Arguments:
    #   index - an index or price S4 'timeSeries' object
    #   spar - a numeric between 0 and, the degree of smoothness
    #   main - main plot title
    #   trace - a logical, should the results be traced?
    #   plot -  a logical, should the results be plotted
    #   at - generate pretty axis positions
    #   format -  a string describing the label format
    
    # Description:
    #   Retroactive Garch11 Volatility Analytics
    
    # FUNCTION:
    
    # Load Library:
    require(fGarch)
    
    # Settings:
    stopifnot(isUnivariate(index))
    if (is.null(main)) main <- "Retroactive Garch Analytics"
    
    # Fit Garch11 Model:
    fit <- fGarch::garchFit(data = 100 * returns(index), trace=trace)
    xseries <- as.timeSeries(fit@data)/100
    Index <- cumulated(xseries)
    colnames(Index) <- colnames(index)
    
    # Turning Points:
    tps <- turnsAnalytics(index=index, spar=spar, trace=trace, doplot=FALSE)
    ablines <- tps$ablines
    
    # Plot Return Series and Standard Deviations:
    xcsd <- timeSeries(data = fit@sigma.t/100, charvec = time(xseries))
    if (doplot) {
      
      # Plot:
      sdPlus <- mean(xseries) + 2 * xcsd
      sdMinus <- mean(xseries) - 2 * xcsd
      range <- range(xseries, sdPlus, sdMinus)
      plot(xseries, 
           main=main, xlab="", ylab=paste("Volatility", colnames(index)), 
           at=at, format=format, 
           # type="l", 
           col="steelblue", ylim=range)
      abline(v = ablines, lty = 3, lwd = 2, col = "grey")
      lines(xseries, col = "steelblue")
      lines(sdPlus, col = "red", lwd = 2)
      lines(sdMinus, col = "red", lwd = 2)
      abline(h = 0, col = "grey", lty = 3)
      box(col = "white")
      box(bty = "l")
      
      # Margin Text:
      # mtext("Volatility Band: 2 sd", adj = 0, side = 4, cex = 0.7, 
      #     col = "darkgrey")
      
      # Add Rmwtrics - Do not Remove!
      mtext("Rmetrics", adj=0, side=4, cex=0.7, col="darkgrey")
      
    }
    
    # Return Value:
    invisible(list(
      index = index, 
      residuals = xseries,
      volatility = xcsd, 
      fit = fit))
  }


###############################################################################


riskmetricsAnalytics <- 
  function(index, spar=0.5, lambda=0.9, main=NULL, 
           trace=TRUE, doplot=TRUE, at=pretty(index), format="%m/%y")
  {    
    # A function implemented by Diethelm Wuertz and Tobias Setz
    
    # Description:
    #   Retroactive riskmetrics analytics
    
    # Arguments:
    #   index - an index or price S4 'timeSeries' object
    #   spar - a numeric between 0 and, the degree of smoothness
    #   main - main plot title
    #   trace - a logical, should the results be traced?
    #   plot -  a logical, should the results be plotted
    #   at - generate pretty axis positions
    #   format -  a string describing the label format
    
    # FUNCTION:
    
    # Settings:
    stopifnot(isUnivariate(index))
    if (is.null(main)) main <- "Retroactive RiskMetrics Analytics"
    
    # Series:
    rets <- returns(index)
    
    # Turning Points:
    tps <- turnsAnalytics(index=index, spar=spar, trace=trace, doplot=FALSE)
    ablines <- tps$ablines
    
    # Indicator:
    sigma <- .emaIndicator(abs(rets), lambda)
    sd <- sqrt(.emaIndicator(rets^2, lambda))
    
    # Plot:
    if(doplot) {
      
      # Plot:
      sdPlus <- mean(rets) + 2 * sd
      sdMinus <- mean(rets) - 2 * sd
      range <- range(rets, sdPlus, sdMinus)
      plot(rets, 
           main="", xlab="", ylab=paste("Volatility", colnames(index)),
           at=at, format=format,
           ylim=range) 
      abline(v=ablines, lty=3, lwd=2, col="grey")
      lines(rets, col="steelblue")
      lines(mean(rets) + 2*sd, col="orange", lwd=2)
      lines(mean(rets) - 2*sd, col="orange", lwd=2)
      lines(mean(rets) + 2*sigma, col="red", lwd=2)
      lines(mean(rets) - 2*sigma, col="red", lwd=2)
      abline(h=0, col="grey", lty=3)
      box(col="white")
      box(bty="l")
      
      # Margin Text:
      # mtext(paste("sd/var Volatility Bands: 2 sd  |  lambda", lambda),
      #   adj=0, side=4, cex=0.7, col="darkgrey")  
      
      # Add Rmwtrics - Do not Remove!
      mtext("Rmetrics", adj=0, side=4, cex=0.7, col="darkgrey")
      
    }
    
    # Return Value:
    invisible(list(
      index = index, 
      analytics = sigma))
  }


###############################################################################


bcpAnalytics <-
  function (index, spar=0.5, FUN=returns, method=c("prob", "mean", "var"),
            main=NULL, trace=TRUE, doplot=TRUE, at=pretty(index), format="%m/%y")
  {
    # A function implemented by Diethelm Wuertz and Tobias Setz
    
    # Description:
    #   Retroactive Bayesian Change Points Analytics
    
    # Arguments:
    #   index - an index or price S4 'timeSeries' object
    #   spar - a numeric between 0 and, the degree of smoothness
    #   main - main plot title
    #   trace - a logical, should the results be traced?
    #   plot -  a logical, should the results be plotted
    #   at - generate pretty axis positions
    #   format -  a string describing the label format
    
    # FUNCTION:
    
    # Load Library:
    require(bcp)
    
    # Settings:
    stopifnot(isUnivariate(index))
    if (is.null(main)) main <- "Retroactive Change Points Analytics"
    
    # Turning Points:
    tps <- turnsAnalytics(index=index, spar=spar, trace=trace, doplot=FALSE)
    positions <- tps$positions
    ablines <- tps$ablines
    
    # BCP Analytics:
    fun <- match.fun(FUN)
    series <- fun(index)
    analytics <- bcp::bcp(series)
    
    # Compose Series:
    method <- match.arg(method)
    select <- c(mean="posterior.mean", var="posterior.var", prob="posterior.prob")
    series(series) <- analytics[[select[method]]]
    
    # Compose y-axis Label:
    ylab <- c(mean="Mean", var="Variance", prob="Probability")
    ylab <- paste("Posterior", ylab[method], colnames(index))
    
    # Select:
    if (method == "prob") {
      ylim <- c(0, 1) 
      prob <- timeSeries(data=analytics$posterior.prob, charvec=time(series))
      prob <- na.omit(prob)
    } else {
      ylim <- range(na.omit(series))
      prob <- NA
    }
    
    # Plot:
    if (doplot) {
      
      # Plot:
      plot(series, type="h", ylim=ylim, las=2, col="grey",
           main=main, xlab="", ylab=ylab,
           at=at, format=format)
      abline(v=ablines, lty=3, lwd=2, col="steelblue")
      points(series, pch=19, cex=0.5)       
      
      box(col = "white")
      box(bty = "l")
      
      # MarginText:
      # mtext(paste("Smooth:", spar), adj = 0, side = 4, cex = 0.7, col = "darkgrey")
      
      # Add Rmwtrics - Do not Remove!
      mtext("Rmetrics", adj=0, side=4, cex=0.7, col="darkgrey")
      
    }
    
    # Return Value:
    invisible(list(
      index = index, 
      analytics = analytics,
      prob = prob))
  }


###############################################################################


pcoutAnalytics <-
  function (index, spar=0.5, main=NULL, 
            trace=TRUE, doplot=TRUE, at=pretty(index), format="%m/%y",
            strong=TRUE, k=2, cs=0.25, outbound=0.25) 
  {
    # A function implemented by Diethelm Wuertz and Tobias Setz
    
    # Description:
    #   Retroactive PCA outlier analytics
    
    # Arguments:
    #   index - an index or price S4 'timeSeries' object
    #   spar - a numeric between 0 and, the degree of smoothness
    #   main - main plot title
    #   trace - a logical, should the results be traced?
    #   plot -  a logical, should the results be plotted
    #   at - generate pretty axis positions
    #   format -  a string describing the label format
    
    # FUNCTION:
    
    # Load Library:
    require(mvoutlier)
    
    # Settings:
    stopifnot(isUnivariate(index))
    if (is.null(main)) main <- "Retroactive Outlier Analytics"
    
    # Turning Points:
    tps <- turnsAnalytics(index=index, spar=spar, trace=trace, doplot=FALSE)   
    positions <- tps$positions
    ablines <- tps$ablines
    
    # Series Settings:
    X <- log(index)
    Y <- returns(index)    
    Z <- cbind(X[-1, ], Y)
    u <- Y
    v <- lag(u, k=-k:k, trim=TRUE)
    U <- u[time(v), ]
    
    # Principal Component Outlier Analytics: 
    ans <- mvoutlier::pcout(v, makeplot=FALSE, explvar=0.99, crit.M1=1/3, 
                            crit.c1=2.5, crit.M2=1/4, crit.c2=0.99, cs=cs, 
                            outbound=outbound)
    
    # Plot: 
    colnames(X) <- paste(colnames(X), "X", sep=":")
    colnames(Y) <- paste(colnames(Y), "Y", sep=":")
    Z <- cbind(X[time(v), ], Y[time(v), ])
    Z@recordIDs <- data.frame(wfinal01=ans$wfinal01, wfinal=ans$wfinal, 
                              wloc=ans$wloc, wscat=ans$wscat, x.dist1=ans$x.dist1, 
                              x.dist2=ans$x.dist2)
    rownames(Z@recordIDs) <- rownames(Z)
    wfinal01 <- ans$wfinal01
    
    if (strong) {
      V <- as.timeSeries(Z@recordIDs)[, "wfinal01"]
      W <- lag(V, k=k, trim=TRUE)
      S <- timeSeries(data=rowSums(W))
      Sfinal <- S == 0
      z <- timeSeries(charvec=time(Z), data=ans$wfinal01)
      z[time(S), ] <- as.integer(!Sfinal)
      wfinal01 <- as.numeric(z)
      Z@recordIDs[, "wfinal01"] <- wfinal01
    }
    
    ans$wfinal <- ans$wfinal01 <- ans$wloc <- ans$wscat <- NULL
    ans$x.dist1 <- ans$x.dist2 <- NULL
    
    U <- Z[, 1]
    wfinal01 <- Z@recordIDs$wfinal01
    datapoints <- length(U)
    
    U <- Z[, 2]
    outliers <- length(U) - sum(wfinal01)
    percent <- round(100 * outliers/length(U), 2)   
    weights <- as.timeSeries(Z@recordIDs)[, "wfinal"]
    
    # Analytics:
    Indicator <- 1 - weights
    invWeights <- as.vector(Indicator)
    extreme = sum(invWeights[invWeights > 0.75]) / length(invWeights)
    analytics <- c(
      mean = mean(1-weights), 
      sd = sd(1-weights), 
      extreme75 = extreme, 
      indicator = Indicator)
    
    ylab <- paste("PC Outlier Prob", colnames(index))
    
    # Plot:
    if(doplot) {
      
      # Plot:
      plot(Indicator, type="h", ylim=c(0, 1), las=2, col="grey", 
           main=main, xlab="", ylab=ylab, 
           at=at, format=format)
      abline(v=ablines, lty=3, lwd=2, col="steelblue")
      points(Indicator, pch=19, cex=0.5)         
      
      # Add Grid and Box:
      grid(NA, ny=NULL) 
      box(col="white")
      box(bty="l")
      
      # Margin Text:
      # mtext(paste("Smooth:", spar), adj=0, side=4, cex=0.7, col="darkgrey")
      
      # Add Rmetrics - Do not Remove!
      mtext("Rmetrics", adj=0, side=4, cex=0.7, col="darkgrey")
      
    }
    
    # Return Value:
    ans <- invisible(list(
      index = index, 
      pcout = ans, 
      series = Z, 
      analytics = analytics,
      prob = Indicator))
    class(ans) <- c("analytics", "list")
    ans
  }


# -----------------------------------------------------------------------------


addRainbow <- 
  function(analytics, palette=rainbow, a=0.3, b=0.8, K=100)
  {
    # A function implemented by Diethelm Wuertz and Tobias Setz
    
    # Example:
    #    analytics <- pcoutAnalytics(index); addRainbow(analytics)
    #    analytics <- bcpAnalytics(index); addRainbow(analytics)
    
    # FUNCTION:
    
    # Get Probability Indicator:
    Indicator <- analytics$prob
    
    # Check:
    stopifnot(isUnivariate(Indicator))
    stopifnot(min(Indicator) >= 0)
    stopifnot(max(Indicator) <= 1)
    
    # Add Spline Smoothed Indicators:
    k.spar <- seq(a, b, length=K)
    for (k in 1:K) {
      curve <- timeSeries::smoothSpline(Indicator, spar=k.spar[k])[, 2]
      if (k == 1) curve1 <- curve
      if (k == K) diff <- curve - curve1
      lines(curve, lwd=2, col=palette(K)[k]) }
    
    # Add Difference Indicator Line:
    lines(diff + 0.5, lwd=2)
    abline(h=0.5, lwd=2, col="orange")
    
    # Return Value:
    invisible()
  }


###############################################################################


waveletSpectrum <-
  function(index, spar=0.5, main=NULL,
           trace=TRUE, doplot=TRUE, at=pretty(index), format="%m/%y") 
  {
    # A function implemented by Diethelm Wuertz and Tobias Setz
    
    # Description:
    #    Morlet Wavelet Analytics
    
    # Arguments:
    #   index - an index or price S4 'timeSeries' object
    #   spar - a numeric between 0 and, the degree of smoothness
    #   main - main plot title
    #   trace - a logical, should the results be traced?
    #   plot -  a logical, should the results be plotted
    #   at - generate pretty axis positions
    #   format -  a string describing the label format
    
    # FUNCTION:
    
    # Load Library:
    require(dplR)
    
    # Settings:
    stopifnot(isUnivariate(index))
    if (is.null(main)) main <- "Morlet Wavelet Spectrum"
    
    # Index Returns:
    returns <- returns(index)
    Time <- 2:length(index)
    Returns <- as.vector(returns)
    
    ans <- dplR::morlet(y1=returns, x1=Time)
    # Returns a list containing:
    # y          Numeric. The original time series.
    # x          Numeric. The time values.
    # wave       Complex. The wavelet transform.
    # coi        Numeric. The cone of influence.
    # period     Numeric. The period.
    # Scale      Numeric. The scale.
    # Signif     Numeric. The significant values.
    # Power      Numeric. The squared power.
    
    # Turning Points:
    tps <- turnsAnalytics(index=index, spar=spar, trace=trace, doplot=FALSE)
    
    # Wavelet Spectrum:     
    p <- as.vector(ans$Power)
    vec <- c(mean=mean(p), sd=sd(p), skew=skewness(p), kurt=kurtosis(p))
    mat <- c(
      O = base::norm(ans$Power, "O"), # "One Norm" 
      I = base::norm(ans$Power, "I"), # "Inf Norm" 
      F = base::norm(ans$Power, "F"), # "Frobenius" 
      M = base::norm(ans$Power, "M")) # "Max Modulus" 
    if (trace) {
      cat("Stats Measures:\n")
      print(vec)
      print(mat)
    }
    analytics <- c(vec, mat)
    
    # Wavelet Parameters:  
    wave.list <- ans  
    wavelet.levels <- quantile(wave.list$Power, probs=seq(from=0, to=1, by=0.1))
    add.coi <- TRUE 
    add.sig <- TRUE
    crn.lab <- "RWI" 
    key.cols <- rev(heat.colors(length(wavelet.levels)))
    key.lab <- expression(paste("Power"^2))
    nyrs=NULL
    crn.col <- "black"
    crn.lwd <- 1
    crn.ylim <- range(wave.list$y) * 1.1
    
    # Settings:
    y <- wave.list$y
    x <- wave.list$x
    wave <- wave.list$wave
    period <- wave.list$period
    Signif <- wave.list$Signif
    coi <- wave.list$coi
    coi[coi == 0] <- 1e-12
    Power <- wave.list$Power
    siglvl <- wave.list$siglvl
    Signif <- t(matrix(Signif, dim(wave)[2], dim(wave)[1]))
    Signif <- Power/Signif
    period2 <- log(period)/log(2)
    ytick <- unique(trunc(period2))
    ytickv <- 2^(ytick)
    coi2 <- log(coi)/log(2)
    coi2[coi2 < 0] <- 0
    coi2.yy <- c(coi2, rep(max(period2, na.rm=TRUE), length(coi2)))
    coi2.yy[is.na(coi2.yy)] <- coi[2]
    yr.vec.xx <- c(x, rev(x))
    par.orig <- par(c("mar", "las", "mfrow"))
    on.exit(par(par.orig))
    nlevels <- length(wavelet.levels)
    key.labs <- formatC(wavelet.levels, digits=4, format="f")
    asp <- NA
    las <- 1   
    xlim <- range(x, finite=TRUE)
    ylim <- range(period2, finite=TRUE)
    ylim[2] <- ylim[2] * 1.1
    
    # Image Plot:
    plot.new()
    plot.window(xlim, ylim, xaxs="r", yaxs="r")
    
    # DW
    # .Internal(filledcontour()) no longer works on 3.0.
    # .Internal(filledcontour(
    #     as.double(x), as.double(period2), Power, 
    #     as.double(wavelet.levels), col=key.cols))
    
    # Use instead:
    graphics::.filled.contour(
      x = as.double(x), 
      y = as.double(period2), 
      z = Power,
      levels = as.double(wavelet.levels), 
      col = key.cols)
    
    title(main=main, xlab="", ylab=paste("Period", colnames(index)))
    box(col="white")
    box(bty="l")
    
    # Add Contours: 
    contour(x, period2, Signif, levels=1, labels=siglvl, 
            drawlabels=FALSE, axes=FALSE, frame.plot=FALSE, 
            add=TRUE, lwd=2, col="black")
    
    # Add Coin of Influence:
    polygon(yr.vec.xx, coi2.yy, density=c(10, 20), 
            angle=c(-45, 45), col="black")
    
    # Add Axis Labels:
    xtick <- NULL
    for (i in 1:length(time(index)))
      xtick <- c(xtick, which.min(abs(time(index) - at[i])))
    axis(1, at=xtick, labels=format(at, format=format))
    axis(2, at=ytick, labels=ytickv, las=2)   
    
    # Add Scale Legend:
    nCol<- length(key.cols)
    positions <- seq(min(x), max(x), length=nCol+1)
    colLevels <- paste(signif(wavelet.levels,2))
    for (i in 1:nCol) {
      lines(x=c(positions[i], positions[i+1]), 
            y=c(1.03,1.03)*max(period2), lwd=3, col=key.cols[i])
      text(x=(positions[i]+positions[i+1])/2, 1.07*max(period2), 
           colLevels[i], cex=0.6) }
    points(positions, rep(1.03*max(period2), length=nCol+1), 
           pch=19, cex=0.7)
    
    # Add Rmetrics - Do not Remove!
    mtext("Rmetrics", adj=0, side=4, cex=0.7, col="darkgrey")
    
    # Return Value:
    invisible(list(
      index = index, 
      spar = spar, 
      wavelet = ans, 
      analytics = list(Time=x, Period=period2, Power=Power)))
  }     


###############################################################################


parAnalytics <- 
  function()
  {
    # A function implemented by Diethelm Wuertz and Tobias Setz
    
    # Description:
    #    Sets the graph frame for an analytics chart
    
    # FUNCTION:
    
    # Graph Frame:
    par(mfrow = c(2, 1))
    par(mar = c(2, 4, 2, 2) + 0.1)
    par(omi = 0.2*c(1, 0.7, 1, 0.7))
    
    # Return Value:
    invisible()
  }


############################################################################### 


