Sample <- function(n, name=NULL) {
  total.counts <- sum(n)
  ra <- n / total.counts

  s <- list(n=n, ra=ra, total.counts=total.counts, name=name)
  class(s) <- "texmex.sample"
  s
}

SampleFit <- function(sample, trunc=TRUE, verbose=TRUE, ...) {
  # fit the observation to a poilog
  if (verbose) {
    message <- sprintf("fitting sample %s%s", ifelse(is.null(sample$name), "", sample$name), ifelse(trunc, " with truncation", ""))
    print(message)
  }
  res <- poilogMLE(sample$n, trunc=trunc, ...)
  
  # compute the empirical pdf and cdf in a table
  # tabulate n+1 so that the resulting vector has abundance of counts=0 in the first slot
  if (trunc) {
    epdf.table <- tabulate(sample$n[sample$n > 0] + 1)
  } else {
    epdf.table <- tabulate(sample$n + 1)
  }
  epdf.values <- epdf.table / sum(epdf.table)
  ecdf.values <- cumsum(epdf.values)
  
  # compute the theoretical pdf and cdf using the poilog fit
  # start at 0 b/c that's where tabulate n+1 starts
  tpdf.values <- dpoilog(0:max(sample$n), res$par['mu'], res$par['sig'], trunc=trunc)
  F.values <- cumsum(tpdf.values)

  # look up the pdf and cdf for every OTU
  # n=0 means 1st entry in tcdf; n=1 means 2nd, etc.
  tpdf <- tpdf.values[sample$n + 1]
  F <- F.values[sample$n + 1]
  
  # compute rescaled reads
  z <- (log(sample$n) - res$par['mu']) / res$par['sig']
  
  # assign and output
  fit <- list(rest=res, epdf.values=epdf.values, ecdf.values=ecdf.values, tpdf.values=tpdf.values, F.values=F.values,
    tpdf=tpdf, F=F, z=z, name=sample$name)
  class(fit) <- "texmex.fit"
  fit
}

PlotFitPP <- function(f, log=FALSE, npoints=10) {
  # PP plot of a texmex fit

  # extend the empirical cdf with 1's
  display.ecdf <- c(f$ecdf.values, rep(1.0, length(f$F.values) - length(f$ecdf.values)))

  xlab <- ""
  ylab <- ""
  
  if (log) {
    # invert the cdf and take the log in xy
    xlab <- "100 - empirical percent"
    ylab <- "100 - theoretical percent"
    plog <- "xy"
    show <- unique(cbind(100*(1.0 - display.ecdf), 100*(1.0 - f$F.values)))
    decr <- TRUE
  } else {
    xlab <- "empirical percent"
    ylab <- "theoretical percent"
    plog <- ""
    show <- unique(cbind(100*display.ecdf, 100*f$F.values))
    decr <- FALSE
  }
  
  lwd <- 2
  cex <- 0.75
  sorted.show <- show[order(show[,1], decreasing=decr), ]
  plot(sorted.show, xlab=xlab, ylab=ylab, log=plog, type='l', cex=cex, lwd=lwd, main=f$name, bty='n')
  
  points(sorted.show[2:npoints, ], cex=cex, pch="|")
  lines(rbind(c(0,0), c(100,100)), lty=3, lwd=lwd, lend=2)
}

SamplePair <- function(sample0, sample1, name=NULL, trunc0=TRUE, trunc1=TRUE) {
  if (class(sample0) != "texmex.sample" || class(sample1) != "texmex.sample") stop("both objects in SamplePair must be Samples")

  # check that the two samples have the same length
  if (length(sample0$n) != length(sample1$n)) stop ("samples in pair must have same length")

  # fit the samples
  fit0 <- SampleFit(sample0, trunc=trunc0)
  fit1 <- SampleFit(sample1, trunc=trunc1)

  # compute the change in the cdf and zscore for all rows
  dF <- fit1$F - fit0$F
  dz <- fit1$z - fit0$z

  # compute the relative abundance statistics
  dra <- sample1$ra - sample0$ra
  lfc <- log10(sample1$ra) - log10(sample0$ra)
  
  sp <- list(sample0=sample0, sample1=sample1, fit0=fit0, fit1=fit1, dF=dF, name=name, dra=dra, lfc=lfc, dz=dz)
  class(sp) <- "texmex.pair"
  sp
}

PlotPair <- function(pair, log='xy', fit=TRUE, highlight=NULL) {
  # show the before/after counts for this sample
  if (log == '') {
    x <- pair$sample0$n
    y <- pair$sample1$n
    if (fit) m <- lm(y ~ x)
  } else {
    rows <- pair$sample0$n>0 & pair$sample1$n>0
    x <- pair$sample0$n[rows]
    y <- pair$sample1$n[rows]   
    if (fit) m <- lm(log(y) ~ log(x)) 
  }
  
  plot(x, y, xlab=pair$sample0$name, ylab=pair$sample1$name, main=pair$name, log=log, asp=1)
  if (!is.null(highlight)) points(pair$sample0$n[highlight], pair$sample1$n[highlight], col='red')
  if (fit) abline(m, lty='dashed')
}

SampleQuad <- function(control, treatment, name=NULL) {
  # check input classes
  if (class(control) != "texmex.pair" || class(treatment) != "texmex.pair") {
    stop(sprintf("input to SampleQuad must be SamplePair objects; instead got %s and %s", class(control), class(treatment)))
  }

  sq <- list(control=control, treatment=treatment, name=name)
  class(sq) <- "texmex.quad"
  sq
}

PlotQuad <- function(quad, highlight=NULL, dF=FALSE, jitter.amount=0.0) {
  # plot the dz against one another, unless dF is specified
  if (dF) {
    x <- quad$control$dF
    y <- quad$treatment$dF
  } else {
    x <- quad$control$dz
    y <- quad$treatment$dz
  }
  
  if (jitter.amount != 0.0) {
    x <- jitter(x, amount=jitter.amount)
    y <- jitter(y, amount=jitter.amount)
  }
  
  if (dF) {
    plot(x, y, xlab=quad$control$name, ylab=quad$treatment$name, main=quad$name, xlim=c(-1.0,1.0), ylim=c(-1.0,1.0))
  } else {
    plot(x, y, asp=1, xlab=quad$control$name, ylab=quad$treatment$name, main=quad$name)
  }
  abline(0, 1, lty=2)
  abline(h=0, lty=2)
  abline(v=0, lty=2)
  
  # find and add the infinite points at the graph edge
  fp <- finitize.points(x, y)
  points(fp$x, fp$y)

  if (!is.null(highlight)) points(fp$x[highlight], fp$y[highlight], col='red')
}

finitize.points <- function(x, y) {
  # put infinite points at the edges of the current plot
  # be sure to run after 'plot'!
  
  x[is.infinite(x) & x<0] <- par("usr")[1]
  x[is.infinite(x) & x>0] <- par("usr")[2]
  y[is.infinite(y) & y<0] <- par("usr")[3]
  y[is.infinite(y) & y>0] <- par("usr")[4]
  data.frame(x=x, y=y)
}