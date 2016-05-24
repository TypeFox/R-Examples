## ms_funs.r --- General functions for multisite data (both Markus/Rachel's and Jay's)
## Author: Stephen J Eglen
## Copyright: GPL
## Mon 10 Sep 2001

## Functions and variables with "mm" in them are mainly for Markus Meister's
## data; those with "jay" in them are for Jay.  Some functions are suitable
## for both types.

## Some of this code requires code from the tcltk package; this is loaded
## by using the DEPENDS: field in the package description.

.plot.mm.s <- function(s, whichcells=NULL,
                      beg=min(unlist(s$spikes), na.rm=TRUE),
                      end=max(unlist(s$spikes), na.rm=TRUE),
                      label.cells = FALSE,
                      use.names=TRUE,
                      show.bursts = FALSE,
                      main=NULL, ylab='Unit',
                      xlab="Time (s)",
                      for.figure=FALSE,
                      show.episodes,
                      episode.y=-0.01,
                      ...) {
  ## Plot the spikes.
  ##
  ## WHICHCELLS is a list of cell numbers to plot; the default is to
  ## plot all of the cells.  If NA is given instead of a channel
  ## reference, just add some whitespace rather than a spike train.
  ##
  ## BEG and END are the time range for which we want to
  ## plot spikes.  When evaluating maxtime, some cells may have no
  ## spikes; their lists get converted to NA in unlist() so those NA
  ## values need removing.  By default, BEG will be the time of the
  ## first spike, and END will be the time of the last spike.
  ## If SHOW.BURSTS is true, we plots the bursts rather than the spikes.
  ## If LABELS.CELLS is true, we write the index of each spike train
  ## in the y-axis; if USE.NAMES is true, we use the cell name, rather than
  ## the index.
  ## If FOR.FIGURE is true, we make a slightly different outline, which
  ## is useful for making the figures.

  if (length(whichcells)>0 && is.numeric(whichcells[1])) {
    ## leave whichcells alone.
    ;
  } else {
    whichcells = .names.to.indexes(names(s$spikes), whichcells, allow.na=TRUE)
  }


  if (is.null(main)) {
    main <- basename(s$file)
  }
     
  N <- length(whichcells)
  ticpercell <- 1/N; deltay <- ticpercell * 0.8;
  yminadd <- ticpercell

  if (show.bursts)
    spikes <- s$spikes
  else
    spikes <- s$spikes

  
  if (for.figure) {
    plot( c(beg, end), c(0,1), type='n',
         yaxt="n",
         bty="n",
         main="",
         xaxt="n",
         xaxs="i", yaxs="i",
         xlab="", ylab="", ...)
    mtext(main, side=3, adj=0, line=0.5)
    
  } else {
    plot( c(beg, end), c(0,1), type='n', bty='n',
         yaxt="n", main=main,
         xlab=xlab, ylab=ylab, ...)
  }
  
  ymin <- 0

  have.bursts <- ( (length(s$allb) > 0) && show.bursts)
  for (cell in whichcells) {
    ts <- spikes[[cell]]                #get spike times.
    n <- length(ts)                     #number of spikes.

    if(n > 0) {
      ## check we have spikes; e.g. cell could be NA, so we wouldn't
      ## need to do anything except clear a bit of space.
      ys <- numeric(n) + ymin
      
      segments(ts, ys, ts, ys+deltay, lwd=0.2) #draw spikes.

      ## simple test to see if bursts have been defined.
      if (have.bursts) {
        burst.times <- s$allb[[cell]]
        if (!is.na(burst.times[1])) {
          ## we have some vald burst info.
          nbursts <- nrow(burst.times)
          ##ys <- rep(ymin+deltay/2, nbursts)

          ## alternate height of bursts so that we can sep adjacent bursts.
          ys <- rep(ymin+deltay/2, nbursts)
          shimmy <- deltay*0.25
          odd <- (1:nbursts) %% 2 == 1
          ys[odd] <- ys[odd] + shimmy

          start.burst <- ts[burst.times[,"beg"]]
          ## for the end of the burst, -1 is needed since if start spike
          ## is 20, and i=3, last spike in burst is 22 (spikes 20, 21, 22)
          end.burst <- ts[ burst.times[,"beg"] + burst.times[,"len"] -1]
          segments(start.burst, ys,
                   end.burst, ys,
                   col="red", lwd=2)
          text(start.burst, rep(ymin+deltay*1.1, nbursts),
               labels=burst.times[,"len"], col="blue")
        }
      }
    }
    ymin <- ymin + yminadd
  }

  if (label.cells) {
    allys <- seq(from=yminadd/2, by=yminadd, length=N)
    if (use.names) {
      labels <- names(spikes)[whichcells]
    } else {
      labels <- whichcells
    }
    ##mtext(labels, side=2, at=allys, las=1)
    ## Try to draw as axis, but still they overlap...
    axis(2, at=allys, labels=labels, las=1, tick=F)
  }

  if (missing(show.episodes)) {
    ## default.
    show.episodes <- ("episodes" %in% names(s))
  }
  
  if (show.episodes) {
    segments(s$episodes[,"beg"], episode.y, s$episodes[,"end"],
             episode.y, col='purple',xpd=NA)
  }

}

.draw.spikes <- function (t, tmin, tmax,
                         ht=1, spike.ht, label, xscale) {

  ## Compulsory args:
  ## T is the set of spike times.
  ## TMIN, TMAX are the min and max times to show.
  ##
  ## HT is the overall ht of the plot.  SPIKE.HT is then the height
  ## of each spike.  (spike.ht should be less than ht.)
  ## LABEL is an optional label to put at the top of each plot.
  ## If XSCALE is given, it should be a vector (lo, hi) indicating the
  ## the scalebar to add -- this is just reusing the x-axis.  Alternatively
  ## if XSCALE is NULL, no scalebar is drawn.
  
  if (missing(spike.ht))
    spike.ht <- 0.7 * ht                #spike ht should be 90% of total ht.


  y.low <- 0.1                          # min y-value, should be [0,1].

  ## throw out spikes outside range [tmin, tmax]
  reject.spikes <-  (t < tmin) | (t > tmax)
  if (any(reject.spikes))
    t <- t[-reject.spikes]
  else
    cat("no spikes outside [tmin,tmax]\n")
  ## set up the plot region, but don't draw any points.
  plot( c(tmin, tmax), c(0, ht),
       bty="n",                         #switch off border
       xlab="", ylab="",                #no labelling of axes.
       xaxt="n", yaxt="n",              #no automatic axes.
       xlim=c(tmin, tmax), ylim=c(0, ht),
       type="n")


  ## We can manually add x and y tics, just for checks...
  if (missing(xscale))
    ## probably don't want xaxis in final version.
    axis(1, at=c(tmin, tmax))             #x-axis
  else {
    ## assume xscale is a 2-d vector providing start and stop time of
    ## scalebar.  tck is the tick length.  labels=FALSE prevents
    ## number labelling of the plot.
    if (!is.null(xscale)) {
      stopifnot(length(xscale)==2)
      axis(1, at=c(xscale[1], xscale[2]), labels=FALSE, tck=0)
    }
  }
  
  ##axis(2, at=c(0, ht))                  #y-axis

  ## for each spike at time t, we draw a line from the point (t,0)
  ## to (t,spike.ht).
  y1 <- y.low + seq(0, by=0, along=t) # vector of zeros, of same length as t.
  y2 <- y1 + spike.ht
  segments(t, y1, t, y2)

  ## optionally label the plot.
  if (!missing(label))
    mtext(label, side=3, adj = 0.02)    #draw label on top axis.
  
}



.summary.mm.s <- function(object, ...) {
  cat(paste("Spike data:", object$file, "\n"))
  cat(paste("NCells", object$NCells, "\n"))
  cat(sprintf("Time (s) [%.3f %.3f]\n", object$rec.time[1], object$rec.time[2]))
}

.read.spikes <- function(reader, ...) {
  ## General-purpose reader around all the xyz.read.spikes() functions.
  ## so that e.g.
  ## .jay.read.spikes(...) can be called as:
  ## .read.spikes(..., reader="jay")
  readers <- c("feller", "iit", "litke", "ncl", "sanger", "sun", "jay", "sql")
  if (reader %in% readers) {
    fn <- paste(reader, ".read.spikes", sep="")
    do.call(fn, list(...))
  } else {
    stop("No such reader:", reader)
  }
}
  
######################################################################
## Jay's functions.
######################################################################



.make.jay.layout <- function(names) {
  ## make the layout for Jay's MEA


  positions <- substring(names, 4,5)
  xlim <- ylim <- c(50, 850)            #now in arrays.R
  spacing <- 100

  cols <- as.integer(substring(positions, 1,1)) * spacing
  rows <- as.integer(substring(positions, 2,2)) * spacing
  pos <- cbind(cols, rows)
  
  rownames(pos) <- names
  array <- 'MCS_8x8_100um'
  
  layout <- list(xlim=xlim, ylim=ylim, spacing=spacing,
                 pos=pos, array=array)

  class(layout) <- "mealayout"

  layout

}



##' Read in the .txt file from Neuroexplorer and create a "spikes" data
##' structure.
##' 
##' Read in .txt file and work out array positions...
##' 
##' 
##' @param filename Name of the text file to be read in.
##' @param ids Optional vector of cell numbers that can be analysed, rather
##' than analysing all electrodes in the recording.  Warning: Not implemented
##' in all readers.
##' @param time.interval Bin width (in seconds) for estimating firing rate.
##' Defaults to 1 second.
##' @param beg Optional start time.
##' @param end Optional end time.
##' @param min.rate Optional minimal firing rate for an electrode to be
##' accepted.
##' @param ... Remaining arguments that are passed to the appropriate reader.
##' @return Return the data structure 's'.
##' @section METHOD: No fancy tricks used here.  If the data file has
##' information about N different spike trains, the file has N (tab-separated)
##' columns.  Each column then gives the time (in seconds?) of each spike.
##' Different columns are of different lengths since typically each cell will
##' have a different number of spikes.
##' 
##' The txt file of spike times can be compressed (with gzip).
##' 
##' read.spikes() is a wrapper around each xyz.read.spikes() function, so that
##' they can all be called just by specifying reader='xyz'.  Current readers
##' are: "feller", "iit", "litke", "ncl", "sanger", "sun", "jay", "sql".
##' 
##' By default, all spikes are read in.  If beg is given, only spikes occuring
##' after this time (in seconds) are kept.  Likewise, if end is given, only
##' spikes occuring before this time (in seconds) are kept.
##' @seealso \code{\link{sanger.read.spikes}}, \code{\link{feller.read.spikes}}
##' @references No references here.
##' @keywords math
##' @examples
##' 
##' data.file <- system.file("examples", "P9_CTRL_MY1_1A.txt",
##'                          package = "IGM.MEA")
##' s <- .jay.read.spikes( data.file)
##' .fourplot(s)
##' s <- .jay.read.spikes(data.file, beg=400, end=700)
##' .fourplot(s)
##' s2 <- .read.spikes(data.file, beg=400, end=700, reader='jay')
##' \dontrun{
##' s <- .jay.read.spikes("~/ms/jay/p9data.txt")
##' .fourplot(s)                             #summary plot.
##' s$mi <- .make.mi(s)
##' .show.prob.t.r(s)                        #conditional distributions.
##' }
##' 
##' \dontrun{.crosscorrplots(s, autocorr=T, tmax=3, nbins=100,
##'                xcorr.nrows=3, xcorr.ncols=3) #plot autocorrs on screen
##' 
##' ## Plotting just one cross-correlogram is a slightly different matter:
##' .xcorr.plot( s$spikes[[1]], s$spikes[[2]], "1 v 2")}
##' 
##' 
##' @export .jay.read.spikes
.jay.read.spikes <- function(filename, ids=NULL,
                            time.interval=1,
                            beg=NULL, end=NULL,
                            min.rate=0) {
  ## Read in Jay's data set.
  ## IDS is an optional vector of cell numbers that should be analysed
  ## -- the other channels are read in but then ignored.

  fp <- gzfile(.file.or.gz(filename), open="r")

  ## todo: in an ideal world, this limit would not be required...
  max.channels <- 200                   #should typically be 64 or less.
  channels <- character(max.channels)
  ## first read in the channel names
  num.channels <- 0
  read.channel.names <- 1
  while(read.channel.names) {
    x<-scan(fp, "", n=1, sep='\t', quiet=TRUE)
    ## If first letter of item is not "c" then assume we have now
    ## reached the timestamps.
    if (tolower(substr(x,1,2)) != "ch") {
      read.channel.names <- 0
      rest <- scan(fp, sep='\t', quiet=TRUE); close(fp)
      ## last element of `rest' is redundant (there is one more TAB that
      ## is not needed at the end of the file), but we need to keep 
      ## x - this is the first element.
      ## File format documented in ~/ms/jay/JAYDATAFORMAT.txt
      times <- c(as.double(x), rest[1:length(rest)-1])
      ntimes <- length(times)
      dim(times) <- c(num.channels, ntimes/num.channels)
      channels <- channels[1:num.channels] #truncate to right size.
      
    } else {
      ## still reading the channel names.
      num.channels <- num.channels + 1
      if (num.channels > max.channels) {
        stop(paste("num.channels has exceeded max.channels"))
      } else {
        channels[num.channels] <- x
      }
    }
  }

  spikes <- apply(times, 1, .jay.filter.for.na)

  if (!is.null(end))
    spikes <- lapply(spikes, .jay.filter.for.max, max=end)

  if (!is.null(beg))
    spikes <- lapply(spikes, .jay.filter.for.min, min=beg)


  
  if (!is.null(ids) ) {
    if (any(ids>length(spikes)))
      stop(paste("some ids not in this data set:",
                 paste(ids[ids>length(spikes)],collapse=" ")))
    
    spikes <- spikes[ids];
    channels <- channels[ids];
  }

  spikes.range <- range(unlist(spikes))
  if (is.null(beg))  beg <-  spikes.range[1]
  if (is.null(end))  end <-  spikes.range[2]
  rec.time <- c(beg, end)
  if (min.rate > 0 ) {
    
    ## Check for inactive channels.
    nspikes <- sapply(spikes,length)
    durn <- diff(rec.time)
    rates <- nspikes/durn
    inactive <- which(rates < min.rate)
    if (any(inactive)) {
      paste("Removing spikes with low firing rates: ",
            paste(inactive, collapse=' '))
      spikes   =   spikes[-inactive]
      channels = channels[-inactive]
    }
    
    
  }

  names(spikes) <- channels
  
  ## Count the number of spikes per channel
  nspikes <- sapply(spikes, length)

  ## meanfiring rate is the number of spikes divided by the (time of
  ## last spike - time of first spike).  
  meanfiringrate <- nspikes/ (end - beg)

  ## Parse the channel names to get the cell positions.

  layout <- .make.jay.layout(channels)
  
  ## temporary test: shuffle electrode positions.
  ## pos <- pos[sample(1:num.channels),]
  
  ## check that the spikes are monotonic.
  .check.spikes.monotonic(spikes)

  rates <- .make.spikes.to.frate(spikes, time.interval=time.interval,
                                beg=beg, end=end)
  
  ## See if we need to shift any units.  this affects only the
  ## visualisation of the units in the movies.  We assume that "shifted"
  ## positions are stored in the file with same name as data file
  ## except that the .txt is replaced with .sps.  Then each line of this
  ## file contains three numbers:
  ## c dx dy
  ## where c is the cell number to move, and dx,dy is the amount (in um)
  ## by which to move the cells.  If you edit the file, this function
  ## must be called again for the new values to be read in.
  ## The shifted positions are used only by the movie functions and
  ## by the function plot.shifted.jay.pos(s) [this shows all units].
  shift.filename <- sub("\\.txt(\\.gz)?$", ".sps", filename)
  unit.offsets <- NULL                  #default value.
  if (file.exists(shift.filename)) {
    updates <- scan(shift.filename)
    ## must be 3 data points per line
    stopifnot(length(updates)%%3 == 0)
    updates <- matrix(updates, ncol=3, byrow=TRUE)
    units <- updates[,1]
    if (any(units> length(spikes))) {
      stop(paste("some units not in recording...",
                 paste(units[units>=length(spikes)],collapse=",")))
    }
    unit.offsets <- layout$pos*0               #initialise all elements to zero.
    unit.offsets[units,] <- updates[,2:3]
  }
  
  
  
  res <- list( channels=channels,
              spikes=spikes, nspikes=nspikes, NCells=length(spikes),
              meanfiringrate=meanfiringrate,
              file=filename,
              layout=layout,
              rates=rates,
              unit.offsets=unit.offsets,
              rec.time=rec.time
               ## TODO: how to return arguments, expanding all vars?
              ##call=match.call()
              )
  class(res) <- "mm.s"

  ## Electrodes are spaced 100um apart in Jay's rectangular array.
  jay.distance.breaks = c(0, 150, 250, 350, 450, 550, 650, 1000)
  res$corr = .corr.index(res, jay.distance.breaks)

  res


}

.jay.filter.for.na <- function(x) {
  ## Truncate the vector X so that trailing NA entries are removed.
  ## This removes the 'empty' spikes at the bottom of each column when
  ## the .txt file is first read in.
  x.na <- which(is.na(x))
  if (any(x.na))
    x[1:x.na[1]-1]
  else
    x
}

.jay.filter.for.max <- function(x, max) {
  ## Any times greater than MAX are removed.
  x.high <- which(x>max)
  if (any(x.high))
    x[1:x.high[1]-1]
  else
    x
}

.jay.filter.for.min <- function(x, min) {
  ## Any times less than MIN are removed.
  ## e.g. .jay.filter.for.min(c(1,2,3,4), 6) should return "nothing?"
  ## .jay.filter.for.min(c(1,2,3,4), 3) returns 3 4
  x.low <- which(x<min)
  if (any(x.low)) 
    x <- x[-x.low]
  x
}

.shuffle.spike.times <- function (s, noise.sd) {
  ## Return new copy of s, with spike trains shuffled to add Gaussian noise
  ## with sd of noise.sd and zero mean.
  spikes <- s$spikes

  add.noise <- function(x, my.sd) {
    ## helper function to add Gaussian noise to a spike train.
    n <- length(x)
    x2 <- sort(x + rnorm(n, sd=my.sd))
  }
  spikes2 <- lapply(spikes, add.noise, noise.sd)
  .check.spikes.monotonic(spikes2)
  s2 <- s                               #make a copy of s
  s2$spikes <- spikes2
  s2
}


.fourplot <- function(s, names=FALSE, raster.time=NULL) {
  ## Simple 2x2 summary plot of an "s" structure.
  ## raster.time should be a 2-vector of (beg, end) time.
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  
  par(mfrow=c(2,2), oma=c(0,0,2,0), las=1)
  
  plot(s$layout, use.names=names)                        #show layout of electrodes.
  .plot.meanfiringrate(s, main='')
  if (is.null(raster.time)) {
    plot(s, main='', label.cells=names, use.names=names)                      #plot the spikes.
  }
  else {
    plot(s, main='', label.cells=names, use.names=names,
         beg=raster.time[1], end=raster.time[2])
  }

  ##   if	(!is.na(s$corr$corr.indexes[1])) {
  if( any(names(s)=="corr")) {
    .plot.corr.index(s, main='')
  }
  mtext(get.project.plate.name(s$file), side=3, outer=T)
}


######################################################################
## Mutual information code.
## taken from Dan Butt's 2002 J Neurosci paper.
.prob.r <- function(s)  {
  ## Given the distance bins, return the probability of finding
  ## two neurons a distance r apart.
  num.distances <- (s$NCells * (s$NCells - 1))/2

  counts <- table(s$dists.bins[which(upper.tri(s$dists.bins))])
  if (sum(counts) != num.distances)
    stop("sum of counts differs from num.distances",
          sum(counts), num.distances)

  ## turn into probability.
  counts/num.distances
}


## Jay and I discussed what "M" should be.  If we have a pair of spike
## trains with 5 spikes in A and 10 spikes in B, should M be 5*10 or
## just the number of dt's which are less than 4 seconds?  (i.e. the
## value of this.m below.)  It should be the latter we think, since
## otherwise if the two spikes are perfectly correlated but very long,
## there will be many dt's that will be greater than 4 seconds.  Dan's
## figure 1b seems to show that each p(t|r) will sum to 1.0 (taking
## into account the bin size; seems to be around 40 bins per 0.5
## second in the plot).

## So, for a pair of spike trains, we compute the cross-correlogram up
## to 4 seconds, and then normalise the histogram so that it sums to
## 1.  This histogram is then binned according to the distance between
## the cell pair.  We then take the average of all histograms for that
## distance bin to compute the overall p(t|r).

.prob.t.cond.r <- function(s, tmax,n.timebins)
{
  ## Return p(t|r) where r is the bin number.
  spikes <- s$spikes
  distance.bins <- s$dists.bins
  n <- s$NCells
  n.distancebins <- length(s$distance.breaks.strings)
  spikepairs <- integer(n.distancebins)
  nhists <- integer(n.distancebins)

  ## Make a matrix to store the histogram for each bin.  Each row
  ## is a histogram for that distance bin.
  allhists <- matrix(0, nrow=n.distancebins, ncol=n.timebins)
  ##dimnames=list("distance bin", "time"))

  hists.rejected <- 0
  ## For each cell pair, compute the histogram of time differences between
  ## spikes, and bin it according to the distance between the cell pair.
  for (a in 1:(n-1)) {
    n.a <- s$nspikes[a]
    for (b in (a+1):n) {
      n.b <- s$nspikes[b]
      bin <- distance.bins[a,b]

      hist <- .hist.ab(spikes[[a]], spikes[[b]], tmax, n.timebins)
      this.m <- sum(hist)

      if (this.m > 0) {
        ## include bin only if there were counts.
        hist <- hist / (this.m)        #normalise by number of comparisons.
        allhists[bin,] <- allhists[bin,] + hist
        nhists[bin] <- nhists[bin] + 1
        spikepairs[bin] <- spikepairs[bin] + this.m
      } else {
        hists.rejected <- 1 + hists.rejected
      }
    }
  }

  if(hists.rejected > 0)
    cat(paste(hists.rejected, "histograms rejected from .prob.t.cond.r\n"))
  
  if ((hists.rejected + sum(nhists)) != (n*(n-1)/2))
    stop(paste("did we compute enough histograms between cell pairs?",
               sum(nhists), (n*(n-1)/2)))

  ## Compute the average histogram for each distance.
  ## Now take the average of each histogram.  This could be done by matrix
  ## multiplication, but this is also simple.
  ## for (i in 1:n.distancebins) {allhists[i,] <- allhists[i,] / nhists[i]}
  allhists <- allhists / nhists

  list(nhists = nhists,
       allhists = allhists,
       spikepairs = spikepairs,
       m=sum(spikepairs),                             # number of counts.
       tmax=tmax,
       n.timebins=n.timebins
       )
}

.prob.t <- function(p.t.cond.r, p.r)  {
  ## Return p(t)
  p.t <- apply(p.t.cond.r, 2, function(x) { drop(x %*% p.r)}) #scalar product
  if (identical(all.equal(sum(p.t),1), FALSE))
    warning("p.t should sum to 1")
  p.t
}



.make.mi <- function(s, tmax=4) {
  ## Return the mutual information.
  ## this includes the bias term,
  ## For Dan's example:
  ## m <- 42000; n.t <-400; n.r <- 9
  ## mi.bias <- ( (n.t * n.r) - n.t - n.r + 1) / ( 2 * m * log(2))
  p.r <- .prob.r(s)
  l <- .prob.t.cond.r(s, tmax, n.timebins=100)
  p.t.cond.r <- l$allhists
  p.t <- .prob.t(p.t.cond.r, p.r)

  if (identical(all.equal(sum(abs(apply(p.t.cond.r, 1, sum))),1),FALSE))
    stop("at least one p(t|r) does not sum to 1")

  (mi <- sum(p.r * apply(p.t.cond.r, 1,
                         function(x) { sum(x * .my.log2(x/p.t)) })) )
  m <- l$m;
  n.t <- length(p.t)
  n.r <- length(p.r)
  mi.bias <- ( (n.t * n.r) - n.t - n.r + 1) / ( 2 * m * log(2))

  mi <- mi - mi.bias                    # subtract bias.

  res <- list(
              mi=mi,
              mi.bias=mi.bias,
              p.r=p.r,
              p.t.cond.r=p.t.cond.r,
              l=l,
              p.t=p.t)
  res
}

.my.log2 <- function(x) {
  ## Take log2(), but change any NaN to 0, since 0log(0) is defined as zero
  ## because x log x -> 0   as x->0.
  res <- log2(x)
  bad <- which(is.infinite(res))
  if (any(bad))
    res[bad] <- 0
  res
}

.show.prob.t.r <- function(s,comment="")  {
  ## Show the p(t|r) distributions.
  ## comment is an optional string to add to the plot.
  ## .make.mi() must have been done first...

  op <- par(no.readonly = TRUE)
  nbins <- length(s$distance.breaks) -1
  if (nbins == 7)
    par(mfrow=c(4,2))                   # jay
  else 
    par(mfrow=c(4,3))                   # MM
  par(mar=c(4,4,2,2))                   #reduce each margin a bit.
  par(oma=c(1,0,0,0))                   #outer margin, 1 line at bottom.

  
  timebin.tmax <- s$mi$l$tmax;
  timebin.n    <- s$mi$l$n.timebins;
  timebin.wid <- timebin.tmax/timebin.n; timebin.min <- 0
  timebin.times <- seq(from=timebin.min+(timebin.wid/2),by=timebin.wid,
                       length=timebin.n)
  
  for (i in 1:nbins) {
    plot(timebin.times, s$mi$p.t.cond.r[i,],
         ##main=paste(s$file,"r bin",i),
         xlab="time (s)",
         ylab=expression(paste("p(", Delta,"t|r)")),
         main=paste(s$distance.breaks.strings[i], "um, n=",s$mi$l$nhists[i]),
         )
    lines( timebin.times[c(1,timebin.n)], c(1/timebin.n, 1/timebin.n),lty=3)
  }
  plot(timebin.times, s$mi$p.t, main="p(t)",
       xlab="time (s)", ylab="p(t)")
  mtext(paste(s$file, date(), "MI",signif(s$mi$mi,4),comment),side=1,outer=TRUE)

  par(op)                               #restore old params.
}

.count.nab <- function(ta, tb, tmax=0.05) {
  ## C routine to count the overlap N_ab (from Wong et al. 1993)
  z <- .C("count_overlap",
          as.double(ta),
          as.integer(length(ta)),
          as.double(tb),
          as.integer(length(tb)),
          as.double(tmax),
          res = integer(1))
  z$res
}



##' Histogram routines to help compute the cross-correlation between a pair of
##' spike trains.
##' 
##' For a pair of spike trains, TA (train a) and TB, these related routines
##' return the count of the number of spikes in B that occur within some time
##' window [-tmax,tmax] of a spike in A.  For .histbi.ab, we return a histogram
##' of dt values from [-tmax,tmax].  For .hist.ab, we ignore the sign of each dt
##' and just return a histogram in the range [0,tmax].  Finally, for .count.nab,
##' we just return the number of dt values found in the range [-tmax,tmax],
##' rather than binning dt into a histogram.
##' 
##' 
##' @aliases .hist.ab .histbi.ab .count.nab .test.count.hist2.nab
##' .test.count.hist.nab .test.hist.ab .test.histograms.versus.r
##' @param ta Vector of spike times, sorted such that lowest is first.
##' @param tb Vector of spike times, sorted such that lowest is first.
##' @param tmax maximum time (in seconds) to bin
##' @param nbins Number of bins in the histogram.  For .histbi.ab, each bin is
##' of width (2*tmax)/nbins.  For .hist.ab, each bin is (tmax)/nbins wide.
##' @return \code{.hist.ab} returns a histogram of counts ignoring sign.
##' \code{.histbi.ab} returns a histogram of counts respecting sign.
##' \code{.count.nab} returns the number of dt values.
##' @section METHOD: For the histogram routines, each bin is of the form [low,
##' high), with the exception of the last bin (for +tmax), which is of the form
##' [tmax-binwid, tmax].  By assuming the spikes are ordered lowest first, the
##' number of spike comparisons is greatly reduced, rather than comparing each
##' spike with A with each spike in B.
##' @seealso Nothing else yet...
##' @references No references here.
##' @keywords math
##' @examples
##' 
##' 
##' stopifnot(isTRUE(all.equal.numeric(
##'   .histbi.ab(c(0), c(-2, -2, 0, 0, 1, 1,1, 1.8,2), tmax=2, nbins=4),
##'   c(2,0,2,5),
##'   check.attributes=FALSE)))
##' stopifnot(identical(TRUE, all.equal.numeric(
##'   .hist.ab(c(0), c(-2, -2, 0, 0, 1, 1,1, 1.8,  2), tmax=2, nbins=4),
##'   c(2,0,3,4),
##'   check.attributes=FALSE)))
##'  
##' .test.hist.ab()
##' 
##' 
## Following examples are not run since they either take a long time
## or require an "s" structure.
##' \dontrun{
##' .test.histograms.versus.r()
##' .test.count.hist.nab()
##' .test.count.hist.nab(s)
##' .test.count.hist2.nab(s)
##' }
##' 
##' @export .hist.ab
.hist.ab <- function(ta, tb, tmax, nbins) {

  ## C routine to bin the overlap time between two spike trains (TA,
  ## TB) into a histogram with NBINS ranging from to TMAX [0,TMAX].
  ## The sign of time differences is ignored.
  z <- .C("bin_overlap",
          as.double(ta),
          as.integer(length(ta)),
          as.double(tb),
          as.integer(length(tb)),
          as.double(tmax),
          res = integer(nbins),
          as.integer(nbins), PACKAGE="IGM.MEA")

  counts <- z$res
  names(counts) <- .hist.make.labels(0, tmax, nbins, right=FALSE)
  counts
  
}

.histbi.ab <- function(ta, tb, tmax, nbins) {
  ## C routine to bin the overlap time between two spikes into a
  ## histogram up to +/- TMAX.  This is a bidirectional version of
  ## .hist.ab, so the sign of time difference matters and the histogram
  ## ranges in [-TMAX,+TMAX]
  
  z <- .C("bin2_overlap",
          as.double(ta),
          as.integer(length(ta)),
          as.double(tb),
          as.integer(length(tb)),
          as.double(tmax),
          res = integer(nbins),
          as.integer(nbins))

  counts <- z$res
  names(counts) <- .hist.make.labels(-tmax, tmax, nbins, right=FALSE)
  counts
}

.hist.make.labels <- function(tmin, tmax, nbins, right=TRUE) {
  ## Make the labels for the histogram bins.
  ## right=TRUE: Each histogram is of the form
  ## (lo, hi], except for the first bin, which is [lo,hi].
  ##
  ## right=FALSE: Each histogram is of the form
  ## [lo, hi), except for the last bin, which is [lo,hi].
  ## This is an internal function that is used from .hist.ab and histbi.ab.
  breaks <- seq(from=tmin, to=tmax, length=nbins+1)
  dig.lab <- 3
  for (dig in dig.lab:12) {
    ch.br <- formatC(breaks, digits = dig, width = 1)
    if (ok <- all(ch.br[-1] != ch.br[-(nbins+1)])) 
      break
  }
  if (right) {
    ## right closed
    labels <- paste("(", ch.br[-(nbins+1)],",", ch.br[-1], "]",sep="")
    substring(labels[1], 1) <- "["
  } else {
    ## left closed
    labels <- paste("[", ch.br[-(nbins+1)],",", ch.br[-1], ")",sep="")
    substring(labels[nbins], nchar(labels[nbins])) <- "]"
  }

  labels
}



.test.histograms.versus.r <- function() {
  ## Test how well my histograms perform against R's binning methods.
  ## Generate some random data points and see how my binning compares to
  ## the binning produced by R's table facility.
  ## If everything goes okay, it should just produce "99" loops.
  ## This is more thorough than the other tests below.
  min.t <- -2.0; max.t <- 2.0; nbins <- 100
  for (i in 1:99) {
    ## Generate some random data.
    r <- rnorm(90000)
    r <- c(r, (numeric(102) + min.t), (numeric(102) + max.t))

    ## In this case, we will assume that all values should fit within the
    ## bounds of the histogram, i.e. we do not need to sort the data.
    
    ## method 1: clip any value outside boundary to boundary value.
    r <- pmin(max.t, pmax(min.t, r))

    ## method 2: reject any value outside boundary.
    ##valid <- which((r >= min.t) & (r <= max.t)); r<- r[valid]

    ## count the number of elements that should be binned.
    count <-  .count.nab(c(0), r,max.t)

    ## t1 - table from R.
    t1 <- table(cut(r, breaks=seq(from=min.t,to=max.t, length=(2*nbins)+1),
                    right=FALSE,include.lowest=TRUE))
    t2 <- .histbi.ab(c(0), r, tmax=max.t, nbins=2*nbins)
    t3 <- .hist.ab(c(0), r,  tmax=max.t, nbins=nbins)
    if(!isTRUE(all.equal(as.numeric(t1),as.numeric(t2)))) {
      print("first test")
      print(t1); print(t2)
      ##error("these are not equal")
    }
      
    if(!all.equal.numeric(sum(t1), count)) {
      print("2nd test")
      print(sum(t1)); print(count)
      ##error("sum and count are unequal")
    }
    bi.cols <- cbind( nbins:1, (nbins+1):(2*nbins))
    bi.sums <- apply(matrix(t2[bi.cols], ncol=2), 1, sum)
    if( !isTRUE(all.equal(as.numeric(t3), bi.sums))) {
      print("third test")
      print(t3); print(bi.sums)
      ##error("t3 and bi.sums are unequal")
    }
    print(i)
  }
  print("all okay")
}
  
.test.hist.ab <- function() {
  ## Test function to check how .hist.ab() compares to R's hist/cut().
  ## If we say spike train A has one spike at time 0, the histogram
  ## produced for comparing spike train A, B will be the same as
  ## binning the spike times of B.
  ## Here we want times in [0,1] to be binned into 4 bins.
  a <- c(0)

  ## either generate random data, or test boundary's explicitly.  Here
  ## the crucial test is whether 0.5 falls in bin 2 or 3.

  ## The histograms produced by my C-code are [low,high), i.e. closed
  ## on the left, open on the right.  To get the same behaviour in cut()
  ## we need to right=F.
  b <- c(0.1,0.2, 0.4, 0.5, 0.9)
  ##b <- runif(100)
  h <- .hist.ab(a, sort(b), 1.0, 4)
  print(h)
  x <- table(cut(b, breaks=c(0,0.25,0.5,0.75,1.0),right=FALSE,
                 include.lowest=TRUE))

  print(x)
  sum(abs(x-h))                         #should be zero.
  
}

.test.count.hist.nab <- function(s) {
  ## For a set of spike trains, check that the C functions
  ## count_overlap and hist_overlap calculate the same values: the
  ## value returned by count_overlap should be the same as the sum of
  ## the histogram returned by hist_overlap.
  spikes <- s$spikes
  
  n <- s$NCells
  dt <- 0.05
  counts <- array(0, dim=c(n,n))
  for ( a in 1:(n-1)) {
    n1 <- s$nspikes[a]
    for (b in (a+1):n) {
      n2 <- s$nspikes[b]
      count <-  .count.nab(spikes[[a]], spikes[[b]],dt)
      counts[a,b] <- count
      this.hist <- .hist.ab(spikes[[a]], spikes[[b]],dt,5)
      if (sum(this.hist) != count) {
        stop(paste("element", a,b, "count", count,
                   "sum", sum(this.hist)))
      }
    }
  }
  ## return the upper triangular array of counts, just in case you want
  ## to examine it.
  counts
}

.test.count.hist2.nab <- function(s) {
  ## For a set of spike trains, check that the C functions
  ## count_overlap and hist_overlap calculate the same values: the
  ## value returned by count_overlap should be the same as the sum of
  ## the histogram returned by hist_overlap.  Furthermore, the
  ## bi_overlap function should be the same.
  
  spikes <- s$spikes
  
  n <- s$NCells
  dt <- 0.05
  nbins <- 10
  counts <- array(0, dim=c(n,n))

  ## which bins go together for the same absolute time delay.
  ## Each column tells you which bins of the bi histogram should be
  ## added to make the one-way histogram.
  bi.cols <- cbind( nbins:1, (nbins+1):(2*nbins))

  for ( a in 1:(n-1)) {
    n1 <- s$nspikes[a]
    for (b in (a+1):n) {
      n2 <- s$nspikes[b]
      count <-  .count.nab(spikes[[a]], spikes[[b]],dt)
      counts[a,b] <- count
      this.hist <- .hist.ab(spikes[[a]], spikes[[b]],dt,nbins)
      ## when doing the bidirectional, must double the number of bins.
      this.hist.bi <- .histbi.ab(spikes[[a]], spikes[[b]],dt,nbins*2)

      ## then work out the sums of the bins that correspond to the
      ## same absolute time differences.  Works only for an even
      ## number of bins.
      bi.sums <- apply(matrix(this.hist.bi[bi.cols], ncol=2), 1, sum)

      ##print(this.hist.bi); print(bi.sums); print(this.hist); stop("stop");

      if (sum(this.hist) != count) {
        stop(paste("element", a,b, "count", count,
                   "sum", sum(this.hist)))
      }
      if (any (this.hist - bi.sums)) {
        print(this.hist)
        print(bi.sums)
        stop(paste("histbi element", a,b, "count", count,
                   "sum", sum(this.hist)))
      }
    }
  }
  ## return the upper triangular array of counts, just in case you want
  ## to examine it.
  ##counts
  NULL
}

.check.similarity <- function(s, tmax=0.001) {
  ## Check to see if two cells have similar spike trains.
  ## Check pair-wise to see the incidence of coincident spiking
  ## (within TMAX seconds of each other).
  ## Return an array,. showing for each cell pair:
  ## i, j, count, n.i, n.j, frac

  ## where i,j are the numbers of the cells; count is the raw count of
  ## coincident spikes; n.i, n.j are the number of spikes in those
  ## trains, and frac is the count/min(n.i, n.j)

  n.cells <- s$NCells
  n.comparisons <- (n.cells * (n.cells-1))/2
  results <- matrix(0, nrow = n.comparisons, ncol = 6) #results array.
  result.line <- 1
  for (i in 1:(n.cells-1)) {
    n.i <- s$nspikes[i]
    for (j in (i+1):n.cells) {
      count <- .count.nab(s$spikes[[i]], s$spikes[[j]], tmax)
      n.j <- s$nspikes[j]
      frac <- count/ min(c(n.i, n.j))
      results[result.line, ] <- c(i, j, count, n.i, n.j, frac)
      result.line <- 1 + result.line
    }
  }
  colnames(results) <- c("i", "j", "count", "n.i", "n.j", "frac")
  results
}

## Global variable that controls whether the xaxis of the xcorr plot
## is shown.
.xcorr.plot.xaxistimes <- FALSE

.xcorr.plot <-  function(spikes.a, spikes.b,
                        plot.label='',
                        xcorr.maxt=4, bi= TRUE,
                        nbins=100,
                        show.poisson=TRUE,
                        autocorr=FALSE, page.label= date(),
                        pause=TRUE,
                        plot=TRUE) {

  ## Produce the cross-correlation of two spike trains, SPIKES.A and SPIKES.B.
  ## PLOT.LABEL is the text to be drawn under the plot.
  ## If BI is true, the histogram is [-XCORR.MAXT, XCORR.MAXT], and we see 
  ## both the negative and positive parts of the correlogram.
  ## page.label is the label to add at the bottom of the page.
  ## To make an autocorrelation, SPIKES.A and SPIKES.B are the same train,
  ## and set AUTOCORR to true.  (For autocorrelation we exclude "self counts",
  ## when a spike is compare to itself.)
  ## If PAUSE is true, during interactive usage, we pause between screenfulls.
  ## (X only, may not work on windows...)
  ## If PLOT is TRUE (default), show the resulting plot.  If FALSE, just return the
  ## cross-correlation values.
  if (bi) {
    x <- .histbi.ab(spikes.a, spikes.b, xcorr.maxt, nbins)
  } else {
    x <-   .hist.ab(spikes.a, spikes.b, xcorr.maxt, nbins)
  }
  
  if (autocorr) {
    ## We are doing auto-correlation, so subtract Ncells from zero bin.
    ## This is to stop counting the time from one spike to itself.
    zero.bin <- floor((nbins/2)+1)
    x[zero.bin] <- x[zero.bin] - length(spikes.a)
    if (x[zero.bin] < 0)
      stop(paste("zero.bin cannot be reduced below zero",
                 x[zero.bin], length(spikes.a)))
  }

  
  dt.wid <- (2*xcorr.maxt)/nbins        #width (in seconds) of one bin.

  ## Normalize to spikes/sec, by dividing the bin count by (Numspikes*bin)
  ## where Numspikes is the number of spikes in the train, and bin is the
  ## time width of each bin.  (From the Neuroexplorer manual.)
  ##  In contrast, if "probability" is required, the normalisation is that
  ## bin counts are dvided by the number of spikes in the spike train.
  x <- x/ (length(spikes.a) * dt.wid)

  max.val <- signif(max(x),2)

  ## Poisson rate is simply the mean firing rate of the other cell.
  ## This calculated as the number of spikes divided by (time of last
  ## spike minus time of first spike.)
  nspikes.b <- length(spikes.b)
  poisson.rate <- nspikes.b/ (spikes.b[nspikes.b] - spikes.b[1])
  

  if (plot) {
    ## Plot the histogram.  type "l" is line, "h" for impulses.
    ## No axes are added here.
    plot(x, ylim=c(0,max.val), type="l",
         bty="n",
         xlab="", ylab="", xaxt="n",yaxt="n")

    ## if we want a y-axis, rather than "max" line...
    want.yaxis <- TRUE
    if (want.yaxis) 
      axis(2, at = c(0, max.val), las=1)
    
    if (show.poisson) {
      lines(c(1, length(x)), c(poisson.rate, poisson.rate), lty=1, col="cyan")
    }
    ## Now annotate the plot with some info.  Plot the info as a central
    ## "tic mark" along the x-axis (which goes from 1 to nbins)
                                        #   axis(1, (nbins/2),
                                        #        labels=c(paste(plot.label,
                                        #          ifelse(want.yaxis, "", paste(" max", max.val)),
                                        #          ##"", signif(poisson.rate,2),
                                        #          sep="")))

    ## put axis at bottom;
    if (.xcorr.plot.xaxistimes) {
      axis(1, c(1, nbins/2, nbins), labels=c(-xcorr.maxt, 0, xcorr.maxt))
    } else {
      axis(1, c(1, nbins/2, nbins), labels=FALSE)
    }

    ## put label at top:
    mtext(plot.label, side=3, cex=par()$cex)

    
    screen.layout <- par()$mfg
    if ( identical(all.equal.numeric(screen.layout[1:2], c(1,1)), TRUE))
      ## Output the page label for only the first plot of the page.
      mtext(page.label, side=1,outer=TRUE)

    if ( identical(all.equal.numeric(screen.layout[1:2],
                                     screen.layout[3:4]), TRUE)
        && ( (names(dev.cur()) == "X11") || (names(dev.cur()) == "windows"))
        && pause)
      ## If we are using a display and the last plot has just been shown,
      ## wait for the user to press RETURN before displaying next page.
      readline("Press return to see next page of plots.")
  } else {
    ## no plotting, just return cross-correlations.
    x
  }

}


.xcorr.restricted <- function(s, a, b,
                             tmin, tmax,
                             plot.label=paste(a,b,sep=":"),
                             show.poisson=TRUE,
                             xcorr.maxt=5, plot=TRUE) {
  ## Compute the cross-correlation just between TMIN and TMAX for two
  ## cells, A and B.  Times are given in seconds.  If TMIN, TMAX
  ## omitted, they default to min,max time respectively.

  if (missing(tmin)) tmin <- min(unlist(s$spikes))
  if (missing(tmax)) tmax <- max(unlist(s$spikes))
  
  ## Instead of plotting, we could just get the result returned to us.
  spikes.a <-s$spikes[[a]]
  spikes.b <-s$spikes[[b]]

  ## remove spikes outside the time range [tmin, tmax]
  rej.a <- which( (spikes.a < tmin) | (spikes.a > tmax))
  if (any(rej.a)) spikes.a <- spikes.a[-rej.a]

  rej.b <- which( (spikes.b < tmin) | (spikes.b > tmax))
  if (any(rej.b)) spikes.b <- spikes.b[-rej.b]

  ## for debugging, just check the range of spikes are as thought.
  ##print(range(spikes.a))
  ##print(range(spikes.b))

  .xcorr.plot(spikes.a, spikes.b,
             xcorr.maxt=xcorr.maxt,
             bi=TRUE, plot.label=plot.label,
             nbins=100,
             autocorr=FALSE, pause=FALSE,
             show.poisson=show.poisson,
             page.label="page label", plot=plot)

}

.crosscorrplots <- function(s, op.file=NULL, tmax=4, nbins=100,
                           autocorr=FALSE,
                           xcorr.ncols=8, xcorr.nrows=14) {
  ## Show all the cross/auto-correlations for the structure.
  ## OP.FILE is the file to output to; If it ends in ".pdf", a PDF is made,
  ## else a postscript file is made.  If OP.FILE is NULL (the default), output
  ## goes to the screen instead, with a pause after each screenfull.
  ## TMAX (defaults to 4 seconds) is the maximum +/-time; NBINS is the
  ## number of bins.
  ## xcorr.nrows and xcorr.ncols controls the dimensions of the array plot.
  
  xcorr.label <- paste(s$file, date(), "tmax [s]", tmax, "nbins", nbins)

  ## If output file ends in ".pdf", make a pdf, else make a postscript file
  if (is.null(op.file)) {
    op <- par(no.readonly = TRUE)
  }
  else {
    if (any(grep ("\\.pdf$", op.file)))
      pdf(file=op.file, width=8.5, height=11)
    else
      postscript(file=op.file)
  }

  par(oma=c(1,0,0,0), mar=c(1.5,1, 0,0)+0.2, tcl=-0.2, mgp=c(0,0,0))
  par(mfrow=c(xcorr.nrows, xcorr.ncols))

  spikes <- s$spikes

  if (autocorr) {
    ncorrs <- s$NCells;
    cell.comparisons <- cbind( 1:ncorrs, 1:ncorrs, 0)
    breaks <- NULL                      #no gaps in the plots
  } else {
    ## more complicated arrangement for cross-correlation
    d <- s$dists;
    d[lower.tri(d, diag=TRUE)] <- NA
    cellpairs <- which(d>=0, arr.ind=TRUE)
    orders <- order(d[cellpairs])
    d2 <- cbind(cellpairs, d[cellpairs])
    ## Now sort them according to smallest distance first.
    cell.comparisons <- d2[orders,]
    ncorrs <- length(orders);
    ## Use the dists.bins matrix to determine when we need to emit a blank
    ## plot.  This helps in viewing the many plots!
    b <- s$dists.bins;  b[lower.tri(b, diag=TRUE)] <- NA
    breaks <- cumsum(table(b))
  }
  for (n in 1:ncorrs) {
    i <- cell.comparisons[n,1]
    j <- cell.comparisons[n,2]
    d <- cell.comparisons[n,3]
    plot.label <- paste(i, ":", j,
                        if (autocorr) {""} else {paste(" d", d)},
                        sep="")

    .xcorr.plot(spikes[[i]], spikes[[j]],
               xcorr.maxt=tmax, bi=TRUE, plot.label=plot.label,
               nbins=nbins,
               autocorr=autocorr,
               page.label=xcorr.label)
    if (any(breaks == n))               #if this is the end of a distance bin
      plot.new()                        #then make a blank plot.
  }
  if (is.null(op.file))
    par(op)
  else
    dev.off()
}

.check.spikes.monotonic <- function(spikes) {
  ## Check to see that all spike times are monotonically increasing.
  ## The counting and histogram routines assumes that spike times
  ## are sorted, earliest spikes first.
  ## .check.spikes.monotonic( list(c(1,3,5), c(1,5,4)))
  results <- sapply( spikes, function(x) { any(diff(x) <0)})
  if (any(results)) {
    stop(paste("Spikes are not ordered in increasing time",
               paste(which(results),collapse=" "),"\n"))
  }
}
  

.mm.spikes.to.bursts <- function(spikes, burst.sep=2) {
  ## Convert spikes to bursts.
  ## burst.sep is the threshold time between spikes for finding bursts.
  ## .spikes.to.bursts(c(1,2,3, 7,8, 11,12,13,14, 19,20, 23,24))
  ## Note: this is too simplistic, and probably not applicable to animals
  ## at older ages where the spike firing could be almost continuous.
  f <- which( diff(spikes) > burst.sep) +1
  spikes[c(1,f)]
}

.list.to.data.frame <- function(l) {
  ## Convert a list of sublists to a data frame.  Each sublist is assumed
  ## to have the same names; each name forms a column.
  ## .list.to.data.frame( list (list(a=3, name="cat", legs=4),
  ##                           list(a=5, name="human", legs=2),
  ##                           list(a=5, name="snake", legs=0)) )
  ## Performs minimal sanity checking.
  ## Master version is in ~/langs/R/list_to_dataframe.R -- edit that version
  ## and then copy back here!!!
  
  if(length(unique(sapply(l, length))) > 1)
    stop("not all list elements are of the same length")
  
  ## Check that the names of each sublist are identical.
  num.sublists <- length(l)
  if (num.sublists > 1) {
    for(i in 2:num.sublists) {
      ## check that sublist 2,3... has same names as sublist 1.
      if (any(!(names(l[[i]]) == names(l[[1]])))) {
        print(names(l[[1]]))
        print(names(l[[i]]))
        stop(paste("different names in sublists", 1,  i))
      }
    }
  }
  names <- names(l[[1]])
  new.list <- list()
  for (num in 1:length(l[[1]])) {
    t <- sapply(l, function(x) {x[[num]]})
    column <- list(t); names(column)[1] <- names[num]
    new.list[[num]] <- column
  }
  d <- data.frame(new.list)
  d
}



######################################################################
# movie-related functions.

.make.animated.gif <- function (x, beg=1,
                               end=dim(x$rates$rates)[1],
                               delay=10,
                               output="anim.gif",
                               delete.frames=TRUE) {
  ## THIS FUNCTION IS NOW DEPRECATED -- USE MAKE.MOVIEFRAMES INSTEAD.
  ##
  ## Loop over each frame, making a temporary .pbm (black/white) and
  ## then convert it to a GIF.  Temporary gif file names are written
  ## as /tmp/ms.movNNNNN.gif where NNNNN is the frame number.  The
  ## frame number normally has leading zeros (e.g. 00050 rather than
  ## 50) so that the frames are ordered correctly by the * wildcard
  ## when creating the animated gif.
  ##
  ## DELAY is the delay (an integer) in 100ths of a second.
  ## WARNING: this works only on Linux, as it requires a ocuple of
  ## external unix programs.!

  stopifnot(.Platform$OS.type=="unix")
  for (i in beg:end) {
    .plot.rate.mslayout(x, i)
    file <- paste("/tmp/ms.mov", formatC(i,width=5,flag="0"), ".gif", sep='')
    dev2bitmap(file="/tmp/ms.mov.pbm", type="pbmraw")
    system(paste("ppmtogif /tmp/ms.mov.pbm >",file, sep=''), ignore.stderr=TRUE)
    ##file <- paste("/tmp/ms.mov", formatC(i,width=5,flag="0"), ".pbm", sep='')
    ##dev2bitmap(file=file, type="pbmraw")
  }

  ## now make the animated gif.
  system(paste("gifsicle --delay=",delay," --loop /tmp/ms.mov*.gif > ",
               output, sep=''))

  ## Have the option to keep or delete the individual frames after
  ## making the movie.

  if (delete.frames)
    system("rm -f /tmp/ms.mov*.gif /tmp/ms.mov.pbm")

}

##' Generate a movie of MEA activity.
##'
##' The mean firing rate of each unit is computed and represented as a circle with
##' the area proportional to the firing rate.  The sequence of frames are then
##' coerced into a movie.
##' 
##' @param x The "s" object.
##' @param beg start time of the movie
##' @param end end time of the movie
##' @param outputdir directory to store the frames (no slash at end).
##' If directory does not exist, it is created first.
##' @param prefix prefix file name for frames
##' @param show.frames Boolean -- do we show the frames on screen as well?
##' @param seconds Boolean: draw the time above the plot?
##' @param delete.first Boolean: delete the outputdir before starting?
##' @param clean.after Boolean: delete the outputdir after finishing?
##' @param anim.delay time (in seconds) delay between frames.  If delay is zero,
##' do not convert movie.
##' 
##' @return NULL.
##' 
##' @author Stephen Eglen
.make.movieframes <- function (x, beg=1,
                              end=dim(x$rates$rates)[1],
                              outputdir=dirname(tempfile()),
                              prefix="mea",
                              show.frames = interactive(),
                              seconds=TRUE,
                              delete.first=TRUE,
                              clean.after=FALSE,
                              anim.delay=5) {

  ## Loop over each frame, making a PNG (mono) file.
  ## The frame number normally has leading zeros (e.g. 00050 rather than
  ## 50) so that the frames are ordered correctly by the * wildcard.
  ## OUTPUTDIR is the directory where the files are to be stored.  This
  ## should not end in a forward slash (/).  If the directory does not
  ## exist, it is created first.
  ## If DELETE.FIRST is true, we delete all the png files in the output
  ## directory before making any new images.
  ## If SECONDS is true, beg,end are interpreted as time in seconds,
  ## not frames.  These times are then first converted into frame numbers.
  ##
  ## If SHOW.FRAMES is true, we view the frames on the screen as well as
  ## writing them to PNGs.
  ##
  ## Once the frames are made, quicktime on PC can then make a movie of these
  ## frames; or on unix, try: "animate -delay 5 mea*png"
  ##
  ## On unix, we can also use "convert" to make the movies
  ## automatically.  WE can do this by setting ANIM.DELAY to the delay
  ## (in 1/100ths of a second) required between frames.
  
  
  if (substring(outputdir, first=nchar(outputdir))=="/")
    stop(paste("outputdir should not end in slash", outputdir))

  if (!file.exists(outputdir))
    dir.create(outputdir)
  
  if (seconds) {
    ## convert beg, end into frames.
    beg <- .time.to.frame(x$rates$times, beg)
    end <- .time.to.frame(x$rates$times, end)
  }
  
  if (delete.first) {
    ## Delete all movie files before making new set.  Best not to use
    ## unlink as it doesn't accept wildcards on DOS.
    files <- list.files(path=outputdir, full.names=TRUE,
                        pattern=paste(prefix,".*\\.png",sep=''))
    if (length(files)>0)
      file.remove(files)
  }

  ## Show the frames.
  for (i in beg:end) {
    if (show.frames)
      .plot.rate.mslayout(x,i)
    
    file <- paste(outputdir, "/", prefix,
                  formatC(i,width=5,flag="0"),
                  ".png", sep='')
    png(file)
    .plot.rate.mslayout(x, i)
    dev.off()
  }

  cat(paste("Movie frames stored in", outputdir, "\n"))
  
  if (anim.delay > 0) {
    ## We want to make an animation...
    cmd = sprintf("cd %s; convert -delay %d %s*png mea.gif",
      outputdir, anim.delay, prefix)
    ##browser()
    system(cmd)
    cat(sprintf("Output file mea.gif created in %s\n", outputdir))

    if (clean.after) {
      ## remove all the temp files up afterwards
      files <- list.files(path=outputdir, full.names=TRUE,
                          pattern=paste(prefix,".*\\.png",sep=''))
      if (length(files)>0)
        file.remove(files)
    }
  } else {
    cat("Convert these to a movie using Quicktime or on Linux:\n")
    cat("convert -delay 20 *png mea.gif\n")
  }
}

.time.to.frame <- function(times, time) {
  ## Given a vector of TIMES, return the index closest to TIME.
  ## Normally, times will be the vector s$rates$times.
  which.min(abs(times-time))
}

.centre.of.mass <- function(s, beg, end, seconds=TRUE,
                           thresh.num=3, thresh.rate=2) {
  ## Find the centre of mass for a set of spikes.

  ## BEG and END are given in seconds (by default), and converted
  ## into frame numbers here. 
  ## A unit is active if its firing rate is above thresh.rate (Hz).
  ## THRESH.NUM is the minimum number of units that must be active to
  ## draw the Centre of mass.
  ##
  ## We return a list with components:
  ## COM -- a 2-d array giving the centre of mass at each timestep
  ## ACTIVE -- list of units that are active.
  ## METHOD -- the method used to compute CoM.

  first.frame <- 
    if (missing(beg)) 1
    else
      if (seconds)
        .time.to.frame(s$rates$times, beg)
      else
        beg
  
  last.frame <-
    if (missing(end)) length(s$rates$times)
    else
      if (seconds)
        .time.to.frame(s$rates$times, end)
      else
        end

  n.frames <- (last.frame+1 - first.frame)
  com <- array(NA, dim=c(n.frames,2))   #(x,y) coords of COM for each frame.

  rownames(com) <- s$rates$times[first.frame:last.frame]
  colnames(com) <- c("com x", "com y")
  index <- 1
  active <- double(0)                     #vector of units that are active.

  for (f in first.frame:last.frame) {
    above <- which(s$rates$rates[f,] > thresh.rate)
    if (length(above) >= thresh.num) {
      com[index,] <- c( mean(s$layout$pos[above,1]), mean(s$layout$pos[above,2]))
      active <- sort(union(active, above))
    }
    index <- index+1
  }

  res <- list(com=com, active=active, method="thresh")
  class(res) <- "mscom"
  res
}

.centre.of.mass.wt <- function(s, beg, end, seconds=TRUE) {

  ## Find the centre of mass for a set of spikes.
  ## Try a weighting factor version.  All cells included.

  ## BEG and END are given in seconds (by default), and converted
  ## into frame numbers here.
  ## Each unit is weighted by dividing its current firing rate
  ## by the overall firing rate.
  ##
  ## We return a list with two components:
  ## COM -- a 2-d array giving the centre of mass at each timestep
  ## ACTIVE -- list of units that are active.
  ## METHOD -- the method used to compute CoM.
  
  first.frame <- 
    if (missing(beg)) 1
    else
      if (seconds)
        .time.to.frame(s$rates$times, beg)
      else
        beg
  
  last.frame <-
    if (missing(end)) length(s$rates$times)
    else
      if (seconds)
        .time.to.frame(s$rates$times, end)
      else
        end

  n.frames <- (last.frame+1 - first.frame)
  com <- array(NA, dim=c(n.frames,2))   #(x,y) coords of COM for each frame.

  rownames(com) <- s$rates$times[first.frame:last.frame]
  colnames(com) <- c("com x", "com y")
  index <- 1

  for (f in first.frame:last.frame) {
    ## weighting factor of each unit i.
    mass.i <- s$rates$rates[f,] / s$meanfiringrate
    mass <- sum(mass.i)
    com[index, 1] <- sum( mass.i * s$layout$pos[,1]) / mass
    com[index, 2] <- sum( mass.i * s$layout$pos[,2]) / mass
    index <- index+1
  }

  res <- list(com=com, active=NULL, method="wt by mean")
  class(res) <- "mscom"
  res
}

.centre.of.mass.wt2 <- function(s, beg, end, seconds=TRUE,
                               thresh.num=3, thresh.rate=5) {

  ## Find the centre of mass for a set of spikes.
  ## Try a weighting factor version, after we first threshold the units
  ## by the number of cells above a firing rate.
  ## BEG and END are given in seconds (by default), and converted
  ## into frame numbers here.
  ## Each unit is weighted by dividing its current firing rate
  ## by the overall firing rate.
  ##
  ## We return a list with two components:
  ## COM -- a 2-d array giving the centre of mass at each timestep
  ## ACTIVE -- list of units that are active.
  ## METHOD -- the method used to compute CoM.
  
  first.frame <- 
    if (missing(beg)) 1
    else
      if (seconds)
        .time.to.frame(s$rates$times, beg)
      else
        beg
  
  last.frame <-
    if (missing(end)) length(s$rates$times)
    else
      if (seconds)
        .time.to.frame(s$rates$times, end)
      else
        end

  n.frames <- (last.frame+1 - first.frame)
  com <- array(NA, dim=c(n.frames,2))   #(x,y) coords of COM for each frame.

  rownames(com) <- s$rates$times[first.frame:last.frame]
  colnames(com) <- c("com x", "com y")
  index <- 1

  for (f in first.frame:last.frame) {
    above <- which(s$rates$rates[f,] > thresh.rate)
    if (length(above) >= thresh.num) {
      ## weighting factor of each unit i.
      mass.i <- s$rates$rates[f,] / s$meanfiringrate
      mass <- sum(mass.i)
      com[index, 1] <- sum( mass.i * s$layout$pos[,1]) / mass
      com[index, 2] <- sum( mass.i * s$layout$pos[,2]) / mass
    }
    index <- index+1
  }

  res <- list(com=com, active=NULL, method="wt by mean")
  class(res) <- "mscom"
  res
}

.colour.com <- function(com) {
  ## Helper routine to parse the Centre of Mass into consecutive periods
  ## of activity.  This function returns a vector labelling each time-step
  ## of the centre of mass to a wave number.
  nrows <- dim(com)[1]
  colours <- integer(nrows)
  wave <- 0
  in.wave <- FALSE
  for (i in 1:nrows) {
    current <- com[i,1]
    if (in.wave) {
      if (is.na(current)) {
        ## wave has ended
        in.wave <- FALSE
        colour <- NA
      } else {
        ## still within the wave
        colour <- wave
      }
    } else {
      ## see if we are now in a wave.
      if (is.na(current)) {
        colour <- NA
      } else {
        ## new wave has started.
        wave <- wave + 1
        colour <- wave
        in.wave <- TRUE
      }
    }
    colours[i] <- colour
  }

  colours
}

.plot.mscom <- function(x, s, colour=TRUE, show.title=TRUE,
                       label.cells=NULL,
                       ##border=FALSE,
                       max.cols=8,
                       rel.cex=1, ...) {
  ## Plot the centre-of-mass using COLOUR if TRUE.
  ## S is optional, but if given, we get to see electrode positions
  ## and the name of the file.
  par.pty <- par()$pty
  par(pty="s")                          #use square plotting region.
  
  if (colour) {
    colours <- .colour.com(x$com)
    nwaves <- max(colours, na.rm=TRUE)
    
    ## Break up the COM into "waves", consecutive times when we have
    ## the Centre of Mass. We then loop over these to plot with a
    ## different colour for each wave.
    ##com.lines <- sapply(1:nwaves, function(i) {x$com[which(colours==i),1:2]})
    com.lines <- lapply(1:nwaves, function(i) {
      valid <- which(colours==i)
      v <- x$com[valid,1:2]
      ##matrix(data=v, ncol=2)
    })
    col.num <- 0;
    first.plot <- TRUE

    for (i in 1:nwaves) {
      c <- com.lines[[i]]
      if (is.null((dim(c))))
        next
      else {
        col.num <- col.num +1;
        if (col.num >= max.cols) col.num <- 1;
      }
      if (first.plot) {
        times <- rownames(x$com)
        title <- paste(ifelse(missing(s), "unknown file", basename(s$file)),
                       x$method,
                       times[1], times[length(times)])
        first.plot <- FALSE
        plot(c, xlab="", ylab="",
             xlim = s$layout$xlim,
             ylim = s$layout$ylim,
             xaxt="n", yaxt="n",
             col=col.num, asp=1, type="l",
             main=ifelse(show.title,title,""))
      } else {
        lines(c, col=col.num)
      }
      ##text(c[1,1], c[1,2], "*", cex=3*rel.cex)
      ## Draw the starting point and add a bit of jitter.
      text(c[1,1]+(20*runif(1)), c[1,2]+(20*runif(1)), "*",
           col=col.num, cex=3*rel.cex)
    }
    ## draw electrode positions if we have them.
    if(!missing(s)) {
      ## if we don't have active list, just draw them all as empty.
      electrode.cols <- rep(0, dim(s$layout$pos)[1]) #default colour of white.
      if (!is.null(x$active))
        electrode.cols[x$active] <- 1   #black for the active ones.

      points(s$layout$pos, pch=21, bg=electrode.cols, cex=rel.cex*0.9,
             lwd=0.4)
    }
    if (!is.null(label.cells)) {
      text(as.numeric(label.cells[,1]),
           as.numeric(label.cells[,2]),
           labels=label.cells[,3], cex=rel.cex*1)
    }

  } else {
    ## let's not bother with colours.
    plot.default(x$com, type="b",asp=1)
  }

  ## restore "pty" parameter.
  par(pty=par.pty)        
}
  

.show.movie <- function(x, beg=1, end,
                       seconds=TRUE,
                       delay=0.03, ...) {
  ## Show a movie within R.
  ## x is the spikes data structure.
  ## BEG is the number of the first frame.
  ## END is the number of the last frame (defaults to the number of
  ## frames to show).
  ## If seconds is true, BEG and END are taken to be time in seconds, rather
  ## than frame numbers.  These are then converted into frame numbers.
  ## delay gives the delay in seconds between frames.
  if (seconds) {
    beg <- .time.to.frame(x$rates$times, beg)
    if (missing(end))
      end <- dim(x$rates$rates)[1]
    else
      end <- .time.to.frame(x$rates$times, end)
  } else {
    if (missing(end))
      end <- dim(x$rates$rates)[1]
  }
       
  for (f in beg:end) {
    .plot.rate.mslayout(x, f, ...)
    Sys.sleep(delay)
  }
}


.make.spikes.to.frate.old <- function(spikes,
                                 time.interval=1, #time bin of 1sec.
                                 frate.min=0,
                                 frate.max=20,
                                 clip=FALSE,
                                 beg=NULL,
                                 end=NULL
                                 ) {

  ## OLD version: gets slow on the Litke data due to large amount of rearranging of
  ## vectors, unlist() etc.
  ##
  ## Convert the spikes for each cell into a firing rate (in Hz)
  ## We count the number of spikes within time bins of duration
  ## time.interval (measured in seconds).
  ##
  ## Currently cannot specify BEG or END as less than the
  ## range of spike times else you get an error from hist().  The
  ## default anyway is to do all the spikes within a data file.

  ## if clips is set to TRUE, firing rate is clipped within the
  ## values frate.min and frate.max.  This is problably not needed.

  spikes.to.rates <- function(spikes, breaks, time.interval) {
    ## helper function.
    h <- hist(spikes, breaks=breaks,plot=FALSE)
    h$counts/time.interval                #convert to firing rate (in Hz)
  }


  spikes.range <- range(unlist(spikes))
  if (is.null(beg))  beg <-  spikes.range[1]
  if (is.null(end))  end <-  spikes.range[2]
  
  time.breaks <- seq(from=beg, to=end, by=time.interval)
  if (time.breaks[length(time.breaks)] < end) {
    ## extra time bin needs adding.
    ## e.g seq(1,6, by = 3) == 1 4, so we need to add 7 ourselves.
    time.breaks <- c(time.breaks,
                     time.breaks[length(time.breaks)]+time.interval)
  }

  rates <- sapply(spikes, spikes.to.rates, breaks=time.breaks,
                   time.interval=time.interval)
  dimnames(rates) <- NULL
  
  ## Now optionally set the upper and lower frame rates if clip is TRUE.
  if (clip)
    rates <- pmin(pmax(rates, frate.min), frate.max)


  ## Do the average computation here.
  ## av.rate == average rate across the array.
  av.rate <- apply(rates, 1, mean)
  ## We can remove the last "time.break" since it does not correspond
  ## to the start of a time frame.
  res <- list(rates=rates,
              times=time.breaks[-length(time.breaks)],
              av.rate=av.rate,
              time.interval=time.interval)
  res
}


.make.spikes.to.frate <- function(spikes,
                                 time.interval=1, #time bin of 1sec.
                                 frate.min=0,
                                 frate.max=20,
                                 clip=FALSE,
                                 beg=NULL,
                                 end=NULL
                                 ) {
  ## Convert the spikes for each cell into a firing rate (in Hz)
  ## We count the number of spikes within time bins of duration
  ## time.interval (measured in seconds).
  ##
  ## Currently cannot specify BEG or END as less than the
  ## range of spike times else you get an error from hist().  The
  ## default anyway is to do all the spikes within a data file.

  ## Note, we need to check for when there are no spikes; this can
  ## happen when examining a subset of spikes, e.g. a well in a multi-well
  ## plate that was not working.
  ## r <- .make.spikes.to.frate(list(), beg=100, end=200, clip=TRUE)
  nspikes <- lapply(spikes, length)
  nelectrodes <- length(nspikes)
  
  ## if clips is set to TRUE, firing rate is clipped within the
  ## values frate.min and frate.max.  This is problably not needed.
  
  spikes.range <- range(unlist(spikes))
  if (is.null(beg))  beg <-  spikes.range[1]
  if (is.null(end))  end <-  spikes.range[2]
  
  time.breaks <- seq(from=beg, to=end, by=time.interval)
  if (time.breaks[length(time.breaks)] < end) {
    ## extra time bin needs adding.
    ## e.g seq(1,6, by = 3) == 1 4, so we need to add 7 ourselves.
    time.breaks <- c(time.breaks,
                     time.breaks[length(time.breaks)]+time.interval)
   }
  nbins <- length(time.breaks) - 1
  
  z <- .C("frate",
          as.double(unlist(spikes)),
          as.integer(nspikes),
          as.integer(nelectrodes),
          as.double(time.breaks[1]), as.double(time.breaks[nbins]),
          as.double(time.interval),
          as.integer(nbins),
          counts = double(nbins*nelectrodes))

  rates <- matrix(z$counts, nrow=nbins, ncol=nelectrodes)

  ## Check if there are any electrodes to process.
  if (nelectrodes > 0) {
    ## Now optionally set the upper and lower frame rates if clip is TRUE.
    if (clip)
      rates <- pmin(pmax(rates, frate.min), frate.max)

    ## Do the average computation here.
    ## av.rate == average rate across the array.
    av.rate <- apply(rates, 1, mean)
  } else {
    av.rate <- rep(NA, nbins)
  }
  ## We can remove the last "time.break" since it does not correspond
  ## to the start of a time frame.
  res <- list(rates=rates,
              times=time.breaks[-length(time.breaks)],
              av.rate=av.rate,
              time.interval=time.interval)
  res
}

.plot.meanfiringrate <- function (s, beg, end, main=NULL, lwd=0.2, ...) {
  ## Plot the mean firing rate over all the cells at each time step.
  ## Can optionally specify the beginning (BEG) and end (END) time, in
  ## seconds.
  
  if (missing(beg)) beg <- s$rates$times[1]
  if (missing(end)) end <- s$rates$times[length(s$rates$times)]

  if (is.null(main))
    main = basename(s$file)
  
  plot(s$rates$times, s$rates$av.rate, type = "h", xlab = "time (s)",
       xlim=c(beg,end), bty="n", lwd=lwd,
       ylab = "mean firing rate (Hz)", main = main, ...)
}

## Simple statistics of the spike trains.

.fano.array <- function(spikes, fano.timebins=c(0.1, 1.0)) {
  ## Compute .fano factor for set of spike trains over a range of
  ## time bins.
  stopifnot(is.list(spikes))
  a <- sapply(fano.timebins, function(t) { .fano.allspikes(spikes, t)})
  rownames(a) <- 1:length(spikes)
  colnames(a) <- fano.timebins
  a
}

.fano.allspikes <- function(spikes, timebin) {
  ## helper function to compute .fano factor of all spike trains
  ## for one time bin.
  sapply(spikes, function(x) {.fano(x, timebin)[4]})
}

.fano <- function(spikes, bin.wid=0.1) {
  ## Compute the .fano factor for one spike train, and for one bin width.

  ## When computing breaks, sometimes the last break comes before the
  ## last spike time, in which case we remove the spikes that come
  ## after the last break.  This should remove only a very small
  ## number of spikes.
  breaks <- seq(from=0, to=ceiling(max(spikes)), by=bin.wid)
  last.break <- breaks[length(breaks)]
  spikes <- spikes[which( spikes <= last.break)]

  h <-  hist(spikes,
             breaks=breaks,
             plot=FALSE, include.lowest=TRUE)

  counts <- h$counts
  counts.mean <- mean(counts); counts.var <- var(counts)
  counts.fano <- counts.var / counts.mean
  res <- c(bin.wid, counts.var, counts.mean, counts.fano)
  names(res) <- c("bin wid", "var", "mean", ".fano")
  res
}

isi <- function(train) {
  ## Compute the ISI for one spike train.
  n <- length(train)
  if (n>1) {
    isi <- diff(train)
  } else {
    isi <- NA                           #cannot compute ISI with 0 or 1 spike.
  }
  isi
}

.cv.isi <- function(train) {
  ## Given a spike train, compute the CV.ISI.
  n = length(train)
  if ( n >1) {
    isi = diff(train)
    isi.mean = mean(isi)
    isi.sd = sd(isi)
    cv = isi.sd / isi.mean
  } else {
    cv = NA                     #cannot compute ISI with 0 or 1 spike.
  }
  cv
}
  

.plot.cumisi <- function(s, xlim=c(0.01, 30)) {

  ## Show the ISI cumulative historgrams.
  ## Each black line is ISI from one channel.
  ## Red line is the mean ISI across the whole array...
  
  show.isi.cdf <- function(spikes, col='black',lwd=1) {
    if (length(spikes)>1) {
      isi1 <- isi(spikes)
      s <- sort(isi1)
      n <- length(s)
      y <- (1:n)/n
      lines(s, y, col=col,lwd=lwd)
    }
  }

  plot(NA, NA, log='x', xlim=xlim, ylim=c(0,1),
       xlab='ISI (s)', ylab='cumulative probability', type='n', bty='n')

  title(basename(s$file))
  res <- sapply(s$spikes, show.isi.cdf)

  ## This is a hacky way to get the average -- since "tt" will be very long...
  tt <- unlist(s$spikes)
  show.isi.cdf(tt, col='red', lwd=2)

  ## Other ideas for plotting the ISI histogram
  ##x = density(isi1)
  ## plot(x, log='x')
  ##x = hist(isi1, breaks=100)
  
}



## store the maximum and minimum firing rate.  Any firing rate bigger
## than this value is set to this value; this prevents the circles
## from overlapping on the plots.  Likewise, anything smaller than the
## minimum is set to the minimum value.
.jay.ms.max.firingrate <- 10
.jay.ms.min.firingrate <- 0.0                  #min firing rate in Hz.

## if electrodes are 100um, each circle can be no bigger than 50um radius,
## else they will overlap.

.jay.ms.max.rad <- 50                    #radius for highest firing rate.
.jay.ms.min.rad <- 2                     #size of smallest rate.


.rates.to.radii <- function(rates) {
  .rates.to.radii.prop.rad(rates)
}

.rates.to.radii.prop.rad <- function(rates) {
  ## Convert the firing rates RATES into radii, such that radius
  ## is proportional to firing rate.

  ## first ensure rates bounded in [min,max]
  rates <- pmax(pmin(rates,.jay.ms.max.firingrate),
                .jay.ms.min.firingrate)
  
  radii <- .jay.ms.min.rad +
    ((.jay.ms.max.rad - .jay.ms.min.rad)* ( rates - .jay.ms.min.firingrate) /
     (.jay.ms.max.firingrate - .jay.ms.min.firingrate))

  radii
}

.rates.to.radii.prop.area <- function(rates) {
  ## Convert the firing rates RATES into radii, such that area of circle
  ## is proportional to firing rate.
  
  ## first ensure rates bounded in [min,max]
  rates <- pmax(pmin(rates,.jay.ms.max.firingrate),
                .jay.ms.min.firingrate)

  min.area <- pi * (.jay.ms.min.rad^2)
  max.area <- pi * (.jay.ms.max.rad^2)
  area <- min.area +
    ((max.area - min.area)* ( rates - .jay.ms.min.firingrate) /
     (.jay.ms.max.firingrate - .jay.ms.min.firingrate))

  radii <- sqrt(area/pi)

  radii
}

## To compare the effect of the two different methods for converting
## rate to radius:
##
##rates <- seq(from=.jay.ms.min.firingrate,to=.jay.ms.max.firingrate, length=100)
##plot(rates, .rates.to.radii.prop.area(rates), type="l",
##     xlab="rate (Hz)", ylab="radius (um)")
##points(rates, .rates.to.radii.prop.rad(rates),pch=19)

## set to TRUE for colour coding of firing rate; FALSE for radius-encoding.
## Colour-coding seems to flicker a lot more than radius coding...
.plot.rate.colour <- FALSE


.plot.rate.mslayout <- function(...) {
  ## Simple wrapper to decide whether to encode firing rate as a radius or
  ## colour of circle.
  if (.plot.rate.colour)
    .plot.rate.mslayout.col(...)
  else
    .plot.rate.mslayout.rad(...)
}

.plot.rate.mslayout.rad <- function(s, frame.num, show.com=FALSE,
                                   show.time=TRUE,
                                   skip.empty=FALSE,
                                   draw.empty=FALSE) {
  ## New version, fixed for Jay's dimensions.
  ## Plot the given frame number in the multisite layout.
  ## The biggest character size is set by jay.ms.max.firingrate.
  ## If SHOW.COM is true, we show the centre of mass as a green dot.
  ## If SKIP.EMPTY is true, any frames where all circles are at min radius
  ## are not drawn.
  ## If SHOW.TIME is true, write the current time above the plot.
  no.small.dots <- FALSE;               #set this to TRUE/FALSE

  radii <-  .rates.to.radii(s$rates$rates[frame.num,])
  
  ## If the radius is zero, R (on unix) still draws a v. small circle
  ## -- is this a bug?  Anyway, replacing zeros (or small values)
  ## with NAs does the trick if you don't want small circles (but they
  ## act as electrode positions which is handy).

  ## extract the unit positions and optionally update them to account
  ## for offsets, so that cells do not overlap on screen.
  xs <- s$layout$pos[,1]; ys <- s$layout$pos[,2]
  if (!is.null(s$unit.offsets)) {
    xs <- xs + s$unit.offsets[,1]
    ys <- ys + s$unit.offsets[,2]
  }

  draw.anything <- TRUE                 #flag - do not change.
  
  if (no.small.dots) {
    min.radius <- .jay.ms.min.firingrate *.jay.ms.max.rad / .jay.ms.max.firingrate
    small.cells <- which(radii < min.radius)
    if (any(small.cells))
      radii[small.cells] <- NA

    if (length(small.cells) == s$NCells) {
      ## nothing to draw.
      draw.anything <- FALSE
    }
      
  }

  if (draw.anything) {
    if (!skip.empty || (any(radii > .jay.ms.min.rad))) {
      symbols(xs, ys,
              fg="black", bg="black",
              circles=radii,
              xaxt="n", yaxt="n", xlab='', ylab='',
              inches=FALSE,
              xlim=s$layout$xlim, ylim=s$layout$ylim,
              main=ifelse(show.time,
                formatC(s$rates$times[frame.num], digits=1, format="f"),
                "")
              )
      if (show.com) {
        com <- .centre.of.mass(s, frame.num, frame.num, seconds=FALSE)
        if(any(com$active))
          ## for retreat, colour them green.
          ##points(com$com, pch=19, col="green")
          ##points(com$com, pch=21, lwd=0.5, bg="grey")
          ## use a triangle for COM, so distinct from open circles
          ## used in other centre of mass plots.
          points(com$com, pch=23, lwd=0.5, bg="white")
      }
    }
  } else {
    ## nothing to draw, so just draw outline.
    if (draw.empty) 
      plot( NA, NA,
           xaxt="n", yaxt="n", xlab='', ylab='',
           xlim=s$layout$xlim, ylim=s$layout$ylim,
           main=formatC(s$rates$times[frame.num], digits=1, format="f"))
  }
}

## Number of colours to have in the firing rate colourmap
## +the colourmap itself.  Reverse the list so that white is low
## and black is high.
.jay.ms.ncols <- 16
.jay.ms.cmap <- rev(gray(0:(.jay.ms.ncols-1)/(.jay.ms.ncols-1)))

.plot.rate.mslayout.col <- function(s, frame.num, show.com=FALSE,
                                   skip.empty=F) {
  ## Colour indicates firing rate.
  ## Plot the given frame number in the multisite layout.
  ## If SHOW.COM is true, we show the centre of mass as a green dot.
  ## If SKIP.EMPTY is true, any frames where all circles are at min radius
  ## are not drawn.

  ## This time, radii are fixed size (e.g. 45um), but colour varies.
  ## radii <- rep(.jay.ms.max.rad, dim(s$layout$pos)[1])
  radii <- rep(45, dim(s$layout$pos)[1])

  cols <- .rates.to.cols(s$rates$rates[frame.num,])

  ## extract the unit positions and optionally update them to account
  ## for offsets, so that cells do not overlap on screen.
  xs <- s$layout$pos[,1]; ys <- s$layout$pos[,2]
  if (!is.null(s$unit.offsets)) {
    xs <- xs + s$unit.offsets[,1]
    ys <- ys + s$unit.offsets[,2]
  }

  symbols(xs, ys,
          fg="black",
          bg=.jay.ms.cmap[cols],
          bty="n",
          circles=radii,
          xaxt="n", yaxt="n", xlab='', ylab='',
          inches=FALSE,
          xlim=s$layout$xlim, ylim=s$layout$ylim,
          main=formatC(s$rates$times[frame.num], digits=1, format="f"))
  if (show.com) {
    com <- .centre.of.mass(s, frame.num, frame.num, seconds=FALSE)
    if(any(com$active))
      points(com$com, pch=19, col="green")
  }
}

.rates.to.cols <- function(rates) {
  bin.wid <- (.jay.ms.max.firingrate - .jay.ms.min.firingrate) /
    .jay.ms.ncols
  
  cols <- floor( (rates-.jay.ms.min.firingrate)/ bin.wid)+1
  ## in case firing rate is outside range of firingrates, limit values
  ## of cols to within 1:jay.ms.ncols.
  cols <- pmax(cols, 1)                 
  cols <- pmin(cols, .jay.ms.ncols)
}

.plot.rate.mslayout.scale <- function(s) {
  ## Draw the scale bar for the plots.

  if (.plot.rate.colour) {
    x <- seq(from=100, to=700, by=100)
    y <- rep(500, length(x))
    rates <- seq(from=.jay.ms.min.firingrate, to=.jay.ms.max.firingrate,
                 length=length(x))
    radii <- rep(45, length(x))
    cols <- .rates.to.cols(rates)
    
    symbols(x, y,
            fg="black",
            bg=.jay.ms.cmap[cols],
            circles=radii,
            xaxt="n", yaxt="n", xlab='', ylab='',
            inches=FALSE,
            xlim=s$layout$xlim, ylim=s$layout$ylim,
            main="legend")
  } else {
    ## show the radius scale bar.
    x <- seq(from=100, to=700, by=100)
    y <- rep(500, length(x))
    rates <- seq(from=.jay.ms.min.firingrate, to=.jay.ms.max.firingrate,
                 length=length(x))
    radii <-  .rates.to.radii(rates)
    
    symbols( x, y,
            fg="black", bg="black",
            circles=radii,
            xaxt="n", yaxt="n", xlab='', ylab='',
            inches=FALSE,
            xlim=s$layout$xlim, ylim=s$layout$ylim,
            main="legend")
  }
  text(x, y-200, labels=signif(rates,digits=2),cex=0.5)

}

.movie.postage <- function(s, tmin, tmax, file="movies.ps") {
  ## Create a postscript file of the firing rate from TMIN to TMAX (given in
  ## seconds).
  postscript(file=file)
  par(mfrow=c(5, 7))
  par(oma=c(0,0,3,0)); par(mar=c(0,1,3,0))
  par(pty="s")                            #produce a square plotting region.
  .show.movie(s, seconds=TRUE, delay=0,
             beg=tmin, end=tmax,
             show.com=TRUE, skip.empty=TRUE)
  .plot.rate.mslayout.scale(s)
  dev.off()
}

.plot.rate.mslayout.old <- function(s, frame.num) {
  ## Plot the given frame number in the multisite layout.
  ## If you want to plot circles rather than disks, change "pch=19"
  ## to "pch=21".  Do `help("points")' for a summary of plot types.
  ## The biggest character size is set by jay.ms.max.firingrate.
  ## xaxt and yaxt control whether or not the axes are plotted.
  
  plot(s$layout$pos[,1], s$layout$pos[,2], pch=19, xaxt="n", yaxt="n",
       cex=pmin(s$rates$rates[frame.num,],.jay.ms.max.firingrate),
       xlab='', ylab='',
       main=formatC(s$rates$times[frame.num], digits=1, format="f"))
}



.op.picture <- function(pos, rates, iteration) {
  ## output a plot of the multisite array activity as a postscript file.
  ps.scale <- 0.5 ### 1.0                      #overall scale factor for plot.
  
  ps.min.x <- 40; ps.min.y <- 40
  ps.wid <-  560 * ps.scale; ps.ht <- 560 * ps.scale;
  ps.max.x <- ps.min.x + ps.wid
  ps.max.y <- ps.min.y + ps.ht
  ps.centre.x <- 0.5 * (ps.min.x + ps.max.x)
  ps.centre.y <- 0.5 * (ps.min.y + ps.max.y)
  

  ps.header <- paste("%!PS-Adobe-3.0 EPSF-3.0\n",
                     "%%Title: Main canvas\n",
                     "%%BoundingBox: ", ps.min.x, " ", ps.min.y, " ",
                     ps.max.x, " ", ps.max.y, "\n",
                     "%%CreationDate: ", date(), "\n",
                     "%%EndComments\n\n",
                     ##"%% /d { 3 1 roll   moveto 10.0 div drawbox} def\n\n",
                     "/d { 3 mul 20 min 0 360 arc fill } def\n\n",
                     "%%EndProlog\n%%Page: 1 1\n",
                     ps.centre.x, " ", ps.centre.y, " translate\n", sep='')

  ps.trailer <- "showpage\n%%Trailer"

  this.rates <- rates[iteration,]
  ncells <- length(this.rates)

  fname <- paste("frame", formatC(iteration, width=4,flag="0"), sep='')
  zz <- file(paste(fname,".ps",sep=''), "w")  # open an output file connection
  cat(ps.header, file = zz)
  for (i in 1:ncells) {
   p <- paste(pos[i,1], pos[i,2], this.rates[i], "d\n")
   cat(p, file = zz)
  }

  cat(ps.trailer, file = zz)
  close(zz)

  system(paste("mypstopnm -pbm ", paste(fname,".ps",sep='')))
  system(paste("ppmtogif ", paste(fname,".pbm",sep=''),">",
               paste(fname,".gif",sep='')))
         
  
  fname
}


.plot.mealayout <- function(x, use.names=FALSE, ...) {
  ## Decide which function to plot the array layout based on number of cells.
  if (nrow(x$pos) < 100) {
    .plot.mealayout.1(x, use.names, ...)
  } else {
    .plot.mealayout.hi(x, use.names, ...)
  }
}
    
  
.plot.mealayout.1 <- function(x, use.names=FALSE, ...) {

  ## Plot the MEA layout.

  pos <- x$pos
  plot(NA, asp=1,
       xlim=x$xlim, ylim=x$ylim,
       bty="n",
       xlab="spacing (\u00b5m)", ylab="", type="n")
  if (use.names)
    text(pos[,1], pos[,2], rownames(pos), ...)
  else
    text(pos[,1], pos[,2], ...)
}

.plot.mealayout.hi <- function(x, use.names=FALSE, ...) {
  ## Plot the MEA layout, high density version
  pos <- x$pos
  plot(pos, asp=1,
       xlim=x$xlim, ylim=x$ylim,
       bty="n",
       xlab="spacing (\u00b5m)", ylab="", pch=20)
}

.spikes.to.ragged.csv <- function(spikes, filename='test.csv') {
  ## SPIKES is a list of vectors of spike times.
  ## FILENAME is the name of the CSV file to create.
  ##
  ## Each spike train in SPIKES is first padded with NAs at the end of
  ## each train, so that each train is the same length.  These spikes
  ## are then written to a CSV file, chaning NAs to blank entries.
  ## This file can then be read in by other programs, or viewed in a
  ## spreadsheet.

  ## find longest spike train and then na.pad every other spike train
  ## to that length.
  max.nspikes <- max(sapply(spikes, length))

  na.pad <-function(x,n) {
    ## Pad the end of vector X with NA s.t. it is of length N.
    l <- length(x)
    n.na <- n - l
    c(x, rep(NA, n.na))
  }

  spikes.pad <- lapply(spikes, na.pad,n=max.nspikes)
  d <- data.frame(spikes.pad)
  write.csv(d, filename, na='', row.names=F)
}




.print.mm.s <- function(x) {
  ## Default print method for a SPIKES data structure.
  cat("MEA spikes\n")
  cat(basename(x$file), "\n")
  cat("nchannels ", x$NCells, "\n")
}
