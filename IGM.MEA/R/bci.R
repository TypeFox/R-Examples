## ~/proj/maps/ci/test_bci.R


.bci <- function(s) {
  ## Compute burst correlation index.
  ## Assume that 'allb' has been created, i.e. burst analysis performed.

  stopifnot(!is.null(s$allb))
  
  beg.end.spikes <- function(spikes, b) {
    ## Given burst info and spike train, get beg, end, duration information.
    ## Helper function.
    beg <- b[,"beg"]
    end <- b[,"end"]
    cbind(beg=spikes[beg], end=spikes[end])
  }

  n <- s$NCells
  beg.end.times <- mapply(beg.end.spikes, s$spikes, s$allb)

  ## how long is each cell active?
  burst.durations <- sapply(beg.end.times, function(m) {
    sum( m[,"end"] - m[,"beg"])})

  ## how many bursts were found for each cell?
  nbursts <- sapply(beg.end.times, nrow)


  beg <- unlist(sapply(beg.end.times, function(m) m[,"beg"]))
  end <- unlist(sapply(beg.end.times, function(m) m[,"end"]))

  first.burst <- c(0, cumsum(nbursts)[-n])
  names(first.burst) <- s$channels

  total.nbursts <- length(beg)
  c.index <- 0:(total.nbursts-1)


  ## check that it all looks okay: cbind(c.index, beg, end)

  z <- .C("bci_calc",
          as.integer(n),
          as.double(beg), as.double(end),
          as.integer(nbursts), as.integer(first.burst),
          as.double(burst.durations),
          res = double(n*n), PACKAGE="IGM.MEA")

  
  m <- matrix(z$res, n, n)
  dists = .make.distances(s$layout$pos, rm.lower=FALSE)
  diag(m) <- NaN

  ## keep summary of burst information.
  duration = diff(s$rec.time)
  prob.bursts = burst.durations/ duration
  .burst.info = cbind(nbursts=nbursts, durn=burst.durations, prob=prob.bursts)

  independent.mean = mean(.burst.info[,"prob"])^2
  
  
  res = list(dists = dists,
    bci=m,
    .burst.info = .burst.info, independent.mean=independent.mean)
  class(res) <- "bci"

  res
}



.plot.bci <- function(x, show.fit=TRUE, log.y=FALSE) {
  ## Plot the BCI curve.

  n <- nrow(x$bci)
  leading.diagonal <- seq(from=1, by=n+1, length=n)
  dist <- as.vector(x$dists)
  corr <- as.vector(x$bci)

  dist <- dist[-leading.diagonal]
  corr <- corr[-leading.diagonal]

  xlim <- c(0, max(dist))
  if (log.y) {
    ylim <- c(0.0001,1)
    log="y"
  } else {
    ylim <- c(0,1)
    log=""
  }
  plot(dist, corr, pch=20, log=log,
     xlab='intercell distance (um)', ylab='BCI', bty='n', 
       las=1, xlim=xlim, ylim=ylim)

  abline(h=x$independent.mean, col='green', lty=2)

  if (show.fit) {
    id <- cbind(dist, corr)
    y.zero <- which(id[,2]==0)
    if (length(y.zero)>0) {
      id <- id[-y.zero,]
      warning(paste("removing", length(y.zero),"zero entries"))
    }
    x <- id[,1]
    y.log <- log(id[,2])
    fit <- lm(y.log ~ x)
    curve(exp(fit$coeff[1]) * exp(x * fit$coeff[2]), add = TRUE,
          from=0)
  }
  
}
