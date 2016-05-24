## --------------------------------------------------------------------------
##
## Plot RVG frequency table
##
## --------------------------------------------------------------------------

plot.rvgt.ftable <- function(x, rows, alpha=0.01, ...)

  ## ------------------------------------------------------------------------
  ## Plot frequencies in table 'x' (histogram).
  ## The plot range is the union of 2*the confidence intervals
  ## and the range of the frequencies.
  ## ------------------------------------------------------------------------
  ## x     : Object of class "rvgt.ftable" containing frequencies
  ## row   : row(s) for which the histogram should be plotted
  ## alpha : draw lines corresponding to confidence intervals
  ## ...   : further graphical parameters
  ## ------------------------------------------------------------------------
{
  ## --- check arguments ----------------------------------------------------

  if (alpha <= 0 || alpha > 0.11)
    stop ("Invalid argument 'alpha'.")

  ## --- read data ----------------------------------------------------------

  ## get table
  table <- x$count
  
  ## number of bins
  m <- ncol(table)
  
  ## get list of rows
  if (missing(rows)) {
    ## use all
    rows <- 1:x$rep
  }
  else {
    if (!is.numeric(rows) || !all(rows>=1) || !all(rows<=x$rep))
      stop ("Invalid argument 'rows'.")
  }

  ## total samplesize
  n <- x$n * length(rows)

  ## break points
  ubreaks <- x$ubreaks

  ## expected probabilities
  p0 <- diff(ubreaks)

  ## get requested frequencies
  if (length(rows)==1) {
    count <- table[rows,]
  }
  else {
    count <- colSums(table[rows,])
  }

  ## --- check probabilities ------------------------------------------------

  ## we have problems to plot the histogram when probabilities are too small.
  ## thus we collapse the corresponding bins.
  
  tol <- 1e-4/n
  too.small <- which(p0<tol)
  
  if (length(too.small)>0) {
    cumcount <- cumsum(count)[-too.small]
    cumcount[length(cumcount)] <- n
    
    count <- c(cumcount[1],diff(cumcount))
    p0 <- p0[-too.small]
    ubreaks <- ubreaks[-too.small]
  }

  ## --- parameters for plot ------------------------------------------------

  ## normalized frequencies
  freq <- count / (n * p0)
  
  ## maximum and minimum frequencies
  ## (we remove NA and NaN to avoid confusing error messages)
  fmax <- max(freq, na.rm=TRUE)
  fmin <- min(freq, na.rm=TRUE)

  ## standard deviation of bin densities
  s <- sqrt(p0*(1-p0)/n) / p0

  ## half length of confidence intervals
  ci <- s * qnorm(alpha/2., lower.tail=FALSE)

  ## limits for plot
  yl <- c(min(fmin, 1.-ci, 1.-2*mean(ci)),
          max(fmax, 1.+ci, 1.+2*mean(ci)))
  xl <- c(0,1)

  ## --- create plot --------------------------------------------------------

  ## create plotting aera with labels
  plot(xl,yl,xlim=xl,ylim=yl,type="n",
       xlab="F(x)", ylab="normalized frequencies", ...)

  ## draw histogram
  polygon( rep(ubreaks,each=2), c(yl[1],rep(freq,each=2),yl[1]), col="light blue", lwd=0.1 )

  ## add expected value
  abline(1,0,col="dark green",lwd=2)
  
  ## add confidence intervals
  ## we treat equidistributed in a special manner
  if (isTRUE(all.equal(p0, rep(1/length(p0),length(p0))))) {
    ## equidistributed
    abline(1-ci,0,col="red",lwd=2,lty=2)
    abline(1+ci,0,col="red",lwd=2,lty=2)
  }
  else {
    ## not equidistributed
    ub <- rep(ubreaks,each=2)
    ub <- ub[-1]; ub <- ub[-length(ub)]
    lines(ub, rep(1-ci,each=2), col="red",lwd=2,lty=1)
    lines(ub, rep(1+ci,each=2), col="red",lwd=2,lty=1)
  }

  ## draw histogram lines
  lines(rep(ubreaks,each=2), c(yl[1],rep(freq,each=2),yl[1]), col="dark blue", lwd=2 )
}

## --------------------------------------------------------------------------

print.rvgt.ftable <- function(x, ...) {
  type <- switch(x$dtype,
                 "cont"  = "continuous distribution",
                 "discr" = "discrete distribution",
                 stop("internal error"))

  cat("\nrvgtest - ftable:\n")
  cat("   distrib type   =",type,"\n");
  cat("   sample size    =",x$rep,"*",x$n,"=",x$n*x$rep,"\n");
  cat("   # break points =",length(x$ubreaks),"\n\n")
}

## --------------------------------------------------------------------------
