## networkspikes.R --- identify and analsyse network spikes
## Author: Stephen J Eglen
## Copyright: GPL
## Sun 28 Jan 2007
## Taking ideas from Eytan & Marom J Neurosci (2006).



##ns.T = 0.003                             #bin time for network spikes
##ns.N = 10                               #number of active electrodes.

## 2007-07-27: Code merged in from second version, temp in
## ~/proj/sangermea/test_ns.R 

##' Compute network spikes
##' 
##' Compute the network spikes in an MEA recording, by averaging over all the
##' electrodes in the array.
##' 
##' To see the mean network spikes after they have computed, just look at the
##' mean object.
##' 
##' If you wish to see the individual network spikes, try .show.ns(ns, ...)
##' where the remaining args are passed to the plot function.
##' 
##' @aliases .compute.ns .show.ns
##' @param s MEA data structure
##' @param ns.T Bin width (in seconds) for counting spikes.
##' @param ns.N Threshold number of active electrodes required to make network
##' spike.
##' @param sur How many bins either side of peak to retain when computing the
##' mean network spike (default 100 bins either side).
##' @param whichcells An optional vector of electrode names.
##' @param plot Set to TRUE to plot network spikes.
##' @param ns A network spike data structure, returned by
##' \code{\link{.compute.ns}}
##' @param ... Other plot arguments to pass to \code{\link{.show.ns}}
##' @return A list with the following elements: \item{counts}{vector giving the
##' number of active electrodes in each bin; this can be very long!}
##' \item{ns.N}{The value of ns.N used.} \item{ns.T}{the value of ns.T used.}
##' \item{mean}{The profile of the mean network spike (this is a time series
##' object)} \item{measures}{If N network spikes were found, this is a matrix
##' with N rows, one per network spike.} \item{brief}{A short vector
##' summarizing the network spikes.  n: number of spikes; peak.m, peak.sd: mean
##' and sd of the peak height; durn.m, durn.sd: mean and sd of the duration of
##' the network spike.}
##' @author Stephen Eglen
##' @seealso \code{\link{sanger.read.spikes}}
##' @references Eytan and Marom (2006) J Neuroscience.
##' @keywords Network spikes, MEA analysis
##' @examples
##' 
##' data.file <- system.file("examples", "TC89_DIV15_A.nexTimestamps",
##'                          package = "IGM.MEA")
##' s <- sanger.read.spikes( data.file, beg=400, end=700)
##' s$ns <- .compute.ns(s, ns.T=0.003, ns.N=10,sur=100)
##' plot(s$ns, ylab='Count', xlab='Time (s)')
##' plot(s$ns, xlim=c(450, 500),
##'      xlab='Time (s)', ylab='Count')
##' 
##' plot(s$ns$mean, xlab='Time (s)', ylab='Count', main='Mean NS')
##' summary(s$ns)
##' s$ns$brief
##' ## .show.ns(s$ns)  # This shows each network spike!  Can take a long time.
##' 
.compute.ns <- function(s, ns.T, ns.N, sur=100, whichcells=NULL,
                       plot=FALSE) {
  ## Main entrance function to compute network spikes.
  ## Typical values:
  ## ns.T: 0.003 (time in seconds)
  ## ns.N: 10
  ## sur: 100

  indexes = .names.to.indexes(names(s$spikes), whichcells, allow.na=TRUE)
  if (length(indexes)==0) {
    ## No electrodes were found matching "whichcells"
    ## so just return brief information summarising no network activity.
    ns <- list()
    ns$brief <- c(n=NA, peak.m=NA, peak.sd=NA, durn.m=NA, durn.sd=NA)
  } else {
    counts <- .spikes.to.count2(s$spikes[indexes], time.interval=ns.T)
    p <- .find.peaks(counts, ns.N)
    ns <- list(counts=counts, ns.N=ns.N, ns.T=ns.T)
    class(ns) <- "ns"
    m <- .mean.ns(ns, p, plot=plot, nrow=4, ncol=4, ask=FALSE, sur=sur)
    if (is.null(m)) {
      ## No network spikes found.
      ns$brief <- c(n=0, peak.m=NA, peak.sd=NA, durn.m=NA, durn.sd=NA)
    } else {
      ns$mean <- m$ns.mean; ns$measures <- m$measures
      peak.val <- ns$measures[,"peak.val"]
      durn <- ns$measures[,"durn"]
      ns$brief <- c(n=nrow(ns$measures),
                    peak.m=mean(peak.val), peak.sd=sd(peak.val),
                    durn.m=mean(durn, na.rm=TRUE), durn.sd=sd(durn, na.rm=TRUE))

    }
  }
  
  ns
}

.spikes.to.count2 <- function(spikes,
                            time.interval=1, #time bin of 1sec.
                            beg=floor(min(unlist(spikes))),
                            end=ceiling(max(unlist(spikes)))
                            )
{
  ## Convert the spikes for each cell into a firing rate (in Hz)
  ## We count the number of spikes within time bins of duration
  ## time.interval (measured in seconds).
  ##
  ## Currently cannot specify BEG or END as less than the
  ## range of spike times else you get an error from hist().  The
  ## default anyway is to do all the spikes within a data file.
  ##
  ## C version, which should replace spikes.to.count
  ## Returns a time series object.

  ## Each bin is of the form [t, t+dt) I believe, as shown by:
  ## .spikes.to.count2(list( c(0, 6.9), c( 2, 4)))
  
  ## time.breaks <- seq(from=beg, to=end, by=time.interval)
  nbins <- ceiling( (end-beg) / time.interval)

  nspikes <- sapply(spikes, length)     #already computed elsewhere!
  
  z <- .C("ns_count_activity",
          as.double(unlist(spikes)),
          as.integer(nspikes),
          as.integer(length(nspikes)),
          as.double(beg), as.double(end), as.double(time.interval),
          as.integer(nbins),
          counts = integer(nbins))
  ## Return counts as a time series.
  res <- ts(data=z$counts, start=beg, deltat=time.interval)

  res
}

IGM.plot.ns <- function(ns, ...) {
  ## Plot function for "ns" class.
  plot(ns$counts, ...)
  abline(h=ns$ns.N, col='red')

  ##peak.times <- times[ ns$peaks[,1]]
  peak.times <- ns$measures[,"time"]
  peak.val   <- ns$measures[,"peak.val"]
  points(peak.times, peak.val, col='blue', pch=19)

}

.summary.ns <- function(ns) {
  ## Summary function for "ns" class.
  n <- ns$brief["n"]
  cat(sprintf("%d network spikes\n", n))
  peak.m <- ns$brief["peak.m"]
  peak.sd <- ns$brief["peak.sd"]


  durn.m <- ns$brief["durn.m"]
  durn.sd <- ns$brief["durn.sd"]
  cat(sprintf("recruitment %.2f +/- %.2f\n", peak.m, peak.sd))
  cat(sprintf("FWHM %.3f +/- %.3f (s)\n", durn.m, durn.sd))
}

.mean.ns <- function(ns, p, sur=100,
                    plot=TRUE, nrow=8, ncol=8, ask=FALSE) {
  ## Compute the mean network spikes, and optionally show the
  ## individual network spikes.

  ## This code does not check to worry if there is a spike right at either
  ## end of the recording.  naughty!

  if (is.null(p)) {
   if (is.null(ns$measures)) {
     # No ns found in this well
#      cat("*** No network spikes found\n")
      return (NULL)
    } else {
      ## use info previously stored in measures.
      p <- ns$measures
    }
  }

  
  if (plot) {
    old.par <- par(mfrow=c(nrow,ncol), mar=c(2.5,1,1,1),ask=ask)
  }
  ave = rep(0, (2*sur)+1)
  npts = length(ns$counts)
  times <- time(ns$counts)
  measures = matrix(NA, nrow=nrow(p), ncol=4)
  colnames(measures) = c("time", "index", "peak.val", "durn")
  n.ns = 0                              #Number of valid network spikes found
  for (i in 1:nrow(p)) {
    peak.i = p[i,"index"]; lo = (peak.i-sur); hi = peak.i+sur

    ## Check that enough data can be found:
    if ( (lo >0) && ( hi < npts) ) {
      n.ns = n.ns + 1

      dat = ns$counts[lo:hi]
      peak.val = dat[sur+1]
      measures[n.ns, "time"] = times[peak.i]
      measures[n.ns, "index"] = peak.i
      measures[n.ns, "peak.val"] = peak.val
      

      if (plot) {
        plot(dat, xaxt='n', yaxt='n', ylim=c(0,60),
             bty='n', type='l',xlab='', ylab='')
        ##abline(v=sur+1)
        max.time <- ns$ns.T * sur
        axis(1, at=c(0,1,2)*sur,
             ##labels=c('-300 ms', '0 ms', '+300 ms'))
             labels=c(-max.time, 0, max.time))

      }
      
      hm = .find.halfmax(dat, peak.n=sur+1, frac=0.5, plot=plot)
      measures[n.ns, "durn"] = hm$durn* ns$ns.T
      if (plot) {
        text <- sprintf("%d durn %.3f",
                        round(peak.val), measures[n.ns, "durn"])
        legend("topleft", text, bty='n')
      }

      ##dat2 = dat;
      ##dat2[1:(hm$xl-1)] = 0;
      ##dat2[(hm$xr+1):((2*sur)+1)] = 0;
      
      ##k = kurtosis(dat2)
      ##measures[n.ns, 1] = k
      ave = ave + dat


    }
  }

  if (n.ns < nrow(p)) {
    ## Some peaks could not be averaged, since they were at either
    ## beg/end of the recording.
    ## So, in this case, truncate the matrix of results to correct
    ## number of rows.
    measures = measures[1:n.ns,,drop=FALSE]
  }
  
  ## now show the average
  if (n.ns > 0) {
    ave = ave/n.ns
    if (plot) {
      plot(ave, xaxt='n', yaxt='n', bty='n', type='l',xlab='', ylab='')
      legend("topleft", paste("m", round(max(ave))), bty='n')
      .find.halfmax(ave)
    }


    ##stripchart(measures[,1], ylab='K', method='jitter', vert=T, pch=19,
    ##main=paste('kurtosis', round(mean(measures[,1]),3)))
    if (plot) {
      stripchart(measures[,"durn"], ylab='durn (s)', method='jitter',
                 vert=TRUE, pch=19,
                 main=paste('FWHM durn', round(mean(measures[,"durn"]),3)))
    }

    if (plot) {
      par(old.par)
    }

  }

  
  ns.mean = ts(ave, start=(-sur*ns$ns.T), deltat=ns$ns.T)

  list(measures=measures, ns.mean=ns.mean)
}

.show.ns <- function(ns, ...) {
  ## Show the individual network spikes after they have been computed.
  ##
  ## This is useful if you don't show the individual network spikes
  ## when they are first iterated over to calculate the mean.
 
  res <- .mean.ns(ns, p=NULL, ...)
  NULL                                  #ignore result
}




.find.peaks <- function(trace, ns.N) {

  ## Peaks are defined as being all elements between two zero entries
  ## (one at start, one at end) in the time series.  An alternate
  ## definiton might be to require some number N of consecutive zero
  ## entries to surround a peak.
   
  max.peaks = 200000

  npts = length(trace)
  
  peaks = matrix(NA, nrow=max.peaks, ncol=2)
  colnames(peaks) <- c("index", "peak.val")
  n = 0

  inside = FALSE;

  for (i in 1:npts) {

    cur = trace[i]    

    if(inside) {
      ## currently in a peak.
      if (cur == 0) {
        ## no longer inside a peak, save results if peak was tall enough.
        inside=FALSE;

        if (peak > ns.N) {
          n=n+1
          if (n > max.peaks) {
            ## oh oh, need more room.
            browser()
          } else {
            peaks[n,] = c(peak.t, peak)
          }
        }
        
      } else {
        ## still in a peak
        if (cur > peak) {
          peak = cur; peak.t = i;
        }
      }
    } else {
      ## currently outside a peak.
      if (cur > 0) {
        inside = TRUE; peak = cur; peak.t = i
      }
    }
  }

  ## tidy up result at end.

  if (n > 0) {
    peaks = peaks[1:n,,drop=FALSE]
  } else {
    ## no peaks found.
    peaks = NULL
  }
}




.find.halfmax.cute <- function(y) {
  ## Given a peak within DAT, find the FWHM.
  ## This is a cute method, but not robust enough -- since it assumes
  ## that the peak is unimodel -- which may not be the case.
  

  x = 1:length(y)
  p.t = 101                             #HACK!!!
  half.max = y[p.t]/2                   #HACK!!!
  f <- approxfun(x, y)
  f2 <- function(x) { f(x) - half.max }
  l <- uniroot(f2, lower=1, upper=p.t)
  r <- uniroot(f2, lower=p.t, upper=length(y))

  segments(l$root, f(l$root), r$root, f(r$root), col='blue')

}



.find.halfmax <- function(y, peak.n=NULL, plot=TRUE, frac=0.5) {

  ## Given a peak somwhere with Y, find the FWHM.
  ##
  ## If PEAK.N is not null, it will be location of the peak -- this is helpful
  ## when there are multiple peaks within one window, and we want to find
  ## the FWHM of the non-major peak.
  ## By default, frac = 0.5, to find the half max.  Change this to some other
  ## value, e.g. 10% to find 10% onset and offset.
  ## 
  ##
  ## This may fail for a few reasons, e.g. not finding half-max values within
  ## the range, running out of data...
  ## all of which should be counted for!

  n = length(y)

  if (is.null(peak.n))
    peak.n = which.max(y)
  
  peak.val = y[peak.n]

  half.max = peak.val * frac
  
  ## Break the data into three segments:

  ## llllllllllllllllllPrrrrrrrrrrrrrrrrr
  ## P is the peak; examine curve to the left (lll) and to the right (rrr) to
  ## find when the peak has decayed to half max.

  left.y = y[1:(peak.n-1)]
  right.y = y[(peak.n+1):n]

  ## When finding the halfmax value in the left and right side, we
  ## have to check that first all of the halfmax value can be found.
  ## e.g. if the peak value is 50 and all values to the left are 45,
  ## there is no value to the left which is under 25, and so the half
  ## max value cannot be computed.
  

  ## Assume the half max point can be found, we interpolate to find
  ## the point, see below.
  
  underhalf.l = which(left.y < half.max)
  if ( any(underhalf.l) ) {
    xl1 = underhalf.l[length(underhalf.l)]   #get last point under halfmax.
    xl2 = xl1+1
    
    yl1 = y[xl1]; yl2 = y[xl2]
    dy = half.max - yl1


    ## see picture.
    ## below, (xl2 - xl1) should equal 1.
    dx = (dy  *(xl2-xl1))/ (yl2-yl1)
    
    xl.half = xl1 + dx
  } else {
    xl.half = NA                        # could not find half-max to left.
  }

  ## Now work on right of curve.  find first point at which y falls below
  ## half max value.
  underhalf.r = which(right.y < half.max)
  if ( any(underhalf.r) ) {
    xr2 = underhalf.r[1] + peak.n
    xr1 = xr2 - 1
    
    yr1 = y[xr1]; yr2 = y[xr2]
    dy = half.max - yr2
    
    dx = dy * (xr1-xr2)/(yr1-yr2)

    ##stopifnot(dx<0)
    xr.half = xr2 + dx
  } else {
    xr.half = NA
  }

  
  if(plot) {
    ##abline(v=xl.half, col='green'); abline(v=xr.half, col='green'); #temp
    abline(h=peak.val * frac, col='red')
    if (! any(is.na(c(xl.half, xr.half)))) {
      ## check first that both half-maxes are valid.
      segments(xl.half, half.max, xr.half, half.max, col='blue')
    }
  }

  list(xl=xl.half, xr=xr.half, durn=xr.half-xl.half)
}

## now interpolate -- hard way
## a <- approx(x, y, n=length(y)*30)
## lines(a)

## amax.x = which.max(a$y)
## points(a$x[amax.x], a$y[amax.x], col='blue', pch=19)

## ## find right side down.
## half.max = max(y)/2
## rx <- which(a$y[-(1:amax.x)]< half.max)[1] + amax.x

## ## find left side down.
## lx <- which(a$y[1:amax.x]< half.max)
## lx <- lx[length(lx)]
## segments(a$x[lx], a$y[lx],  a$x[rx], a$y[rx], col='blue')

## The "R" way of interpolating -- nice!


.check.ns.plot <- function(counts, p, xlim, ns.N) {

  plot(counts$times, counts$sum, type='l', xlim=xlim,
       xlab="time (s)", ylab='num active channels')
  points(counts$times[p[,1]], p[,2], pch=19, col='blue')
  abline(h=ns.N, col='red')               #threshold line.
}

.ns.bin.peak <- function(p, nbins=12, wid=5) {
  ## Bin values in P into a set of NBINS bins, of size WID.
  ## Bins are right-closed (except for first bin, closed at both ends).
  ## Labels are added onto the bins.
  ##
  ## x <- c(0, 4,5, 20, 54,55, 60)
  ## .ns.bin.peak(x, wid=10, nbins=7 )
  ##

  if ( is.null(p) ) {
    ## no valid values, so no need to make the histogram.
    ## This happens when there are no network spikes.
    p <- 0; invalid <- TRUE
  } else {
    invalid <- FALSE
  }
  
  b <- seq(from=0, by=wid, length=nbins+1)
  max.allowed <- max(b)
  if ( any( above <- which(p > max.allowed)) ) {
    stop("some values above max.allowed")
  }
  h <- hist(p, plot=FALSE, breaks=b)
  c <- h$counts

  if (invalid) {
    ## no valid counts, so set all counts to zero.
    c <- c*0
  }
  
  l <- .hist.make.labels(0, max.allowed, nbins)
  names(c) <- l
  c
}



.ns.identity <- function(s, w=0.1) {
  ## Return the "NSID" matrix, Network Spike IDentity.
  ## Which channels contributed to which network spikes?
  ## W is window of spike identity, +/- 0.1s by default.

  ## peak.times here should be the middle of the NS bin.
  peak.times <- s$ns$measures[,"time"] + (s$ns$ns.T/2)

  ## We do the transpose here so that one row is one network spike.
  nsid <- t(.ns.coincident(peak.times, s$spikes, w))

  
}

.ns.coincident <- function(a, bs, w) {
  ## A is a vector of reference times, sorted, lowest first.  (Here
  ## the times of the peaks of the network spikes.)
  ## B is a list of vectors of spike times.  (Each vector of spike
  ## times is sorted, lowest first.)

  ## For each spike train in B, we see if there was a spike within a
  ## time window +/- W of the time of each event in A.  If there was a
  ## "close" spike in spike train j from B to event i , then
  ## MAT[j,i]=1.
  ## (MAT is a matrix of size CxN, where C is the number of spike
  ## trains in BS, and N is the number of events in A.)
  ## MAT is transposed by the higher level function -- it is kept this
  ## way for ease of the C implementation.
  spike.lens <- sapply(bs, length)
  num.channels <- length(spike.lens)
  z <- .C("coincident_arr",
          as.double(a), as.integer(length(a)),
          as.double(unlist(bs)), as.integer(spike.lens),
          as.integer(num.channels),
          close = integer(length(a)*num.channels),
          as.double(w))
  
  mat <- matrix(z$close, nrow=num.channels, byrow=TRUE)
  dimnames(mat) <- list(channel=1:num.channels, ns.peak=a)
  mat

}

######################################################################
## End of functions
######################################################################

