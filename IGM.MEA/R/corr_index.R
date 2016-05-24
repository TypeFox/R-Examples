## Compute the correlation index, as defined by Wong et al (1993).
## Author: Stephen J Eglen
## Copyright: GPL
## Sun 04 Mar 2007

.corr.index <- function(s, distance.breaks,
                       dt=getOption("IGM.MEA.corr.dt", default=0.05),
                       min.rate=0,
                       corr.method = getOption("IGM.MEA.corr.method", default="CI")) {
  ## Make a correlation index object.
  ## MIN.RATE: if greater than zero, we analyse only spike trains whose
  ## firing rate is greater than this minimal rate.
  ## corr.method is which method to use.
  dists = .make.distances(s$layout$pos)
  dists.bins = .bin.distances(dists, distance.breaks)

  spikes = s$spikes
  if (length(spikes) > 1) {
    ## SJE: 2010-03-17 -- try new version of corr index.
    ##corr.indexes = .make.corr.indexes(spikes, dt, min.rate)
    corr.indexes = NULL
    if (corr.method == "CI") {
      corr.indexes = .make.corr.indexes2(spikes, dt, min.rate)
    }
    if (corr.method == "Tiling") {
      corr.indexes = .tiling.allpairwise(s, dt)
    }
    if (is.null(corr.method)) {
      stop("Corr index not calculated")
    }
    
    
    corr.id = cbind(dist=.my.upper(dists), corr=.my.upper(corr.indexes),
      dist.bin=.my.upper(dists.bins))
    ##corr.id.means = .corr.get.means(corr.id)
    dist.mids = diff(distance.breaks)/2 +
      distance.breaks[-(length(distance.breaks))]
    corr.id.means = .corr.get.means(corr.id, dist.mids)
  } else {
    corr.indexes = NA
    corr.id = NA
    corr.id.means = NA
  }

  ## distance.breaks.strings used only by Mutual Information Code?
  distance.breaks.strings =
    levels(cut(0, distance.breaks, right=FALSE, include.lowest=TRUE))

  res = list(
    ##dists=dists, dists.bins = dists.bins,
    ##corr.indexes = corr.indexes,
    dt = dt,
    corr.id = corr.id,
    corr.id.means = corr.id.means,
    distance.breaks=distance.breaks,
    distance.mids=dist.mids,
    distance.breaks.strings=distance.breaks.strings,
    method=corr.method)

  res
}

.make.distances.check <- function() {
  ## Simple function to check our new and old functions are the same.

  ## For each test, generate synthetic data and check for same results.
  for (i  in 1:10) {
    n = 1000
    pos = matrix(runif(n*2, max=500), n, 2)
    print(system.time(d1 <- .make.distances(pos)))
    print(system.time(d2 <- .make.distances.old(pos)))
    print(all.equal(d1,d2))
  }
}

.make.distances <- function(posns, rm.lower=TRUE) {
  ## POSNS should be a (N,2) array.  Returns a NxN upper triangular
  ## array of the distances between all pairs of cells.

  x <- posns[,1]; y <- posns[,2]
  d = round(sqrt(outer(x, x, "-")^2 + outer(y, y, "-")^2))
  if (rm.lower)
    d[lower.tri(d)] = 0
  
  d
}

.make.distances.old <- function(posns) {
  ## TODO: this is slow for large numbers of electrodes!
  ##
  ## POSNS should be a (N,2) array.  Returns a NxN upper triangular
  ## array of the distances between all pairs of cells.

  ## Currently store distances to the nearest micron, so that it makes
  ## the job of binning distances easier when computing the mean of
  ## correlation index for each "distance".  In Figure 9 of the
  ## Meister 1991 paper, distances are binned into 20um bins to get
  ## round this problem.

  n <- dim(posns)[1]
  dists <- array(0, dim=c(n,n))
  for ( a in 1:n-1)
    for (b in (a+1):n) {
      delta <- posns[a,] - posns[b,]
      dists[a,b] <- round(sqrt( sum(delta**2)))
    }

  dists
}

.bin.distances <- function(dists, breaks) {
  ## DISTS is a upper NxN array.
  ## breaks is a vector of breakpoints.
  ## Return an array of the same size where each distance value is
  ## given a corresponding bin number.

  ## e.g.
  ## dists <- matrix( c(0,0,0, 400,0,0, 80, 1000, 0), nrow=3)
  ## jay.bin.distances(dists)
  ## This binning procedure can then be checked by comparing the
  ## distances and their bins:
  ## plot(.my.upper(s$dists.bins), .my.upper(s$dists))
  ## boxplot(.my.upper(s$dists)~ .my.upper(s$dists.bins))
  
  distances <- .my.upper(dists)
  ## These breaks are hardcoded.

  ##data <- c(0, 100, 700, 900, 400)

  ## Make each bin [low, high) with exception that highest bin is
  ## [low,high]. Labels is false so that we just return numeric count
  ## of bin, rather than a factor.
  bins <- cut(distances, breaks, right=FALSE,
              include.lowest=TRUE, labels=FALSE)
  invalid <- is.na(bins)
  if (any(invalid))
    stop(paste("distances not binned:",
                  paste(distances[which(invalid)],collapse=" ")))
  n <- dim(dists)[1]
  res <- matrix(0, nrow=n, ncol=n)
  res[which(upper.tri(res))] <- bins
  
  res
}


.plot.corr.index <- function(s, identify=FALSE,
                            main=NULL,
                            show.method=TRUE,
                            dot.col='red',
                            show.fit=TRUE, show.ci=FALSE,
                            show.pts=NULL,
                            ylabel="correlation",
                            xlabel=expression(paste("intercell distance (",
                                mu, "m)")),
                            ...) {
  ## Plot the correlation indices as a function of distance.
  ## If identify is TRUE, we can locate cell pairs on the plot using
  ## left mouse button.

  ## Use 'log=y' as one of the extra args if the y-axis should be
  ## drawn on a log scale.
  ## DOT.COL: colour of each dot.
  ## If SHOW.FIT is true, draw the expoential fit.
  ## If SHOW.CI is true, draw the confidence intervals estimated every
  ## 100 um or so.
  ## SHOW.PTS: if TRUE, show individual CI values.  If NULL, the value
  ## is assumed TRUE iff number of cells recorded is less than 100.
  
  dists = s$corr$corr.id[,"dist"]
  corrs = s$corr$corr.id[,"corr"]
  
  if (all(is.na(corrs))) {
    ## no correlation data to show, so just up empty plot.
    plot(NA, xlim=range(dists), ylim=c(1,10),
         xlab=xlabel, ylab=ylabel,
         main=paste(get.project.plate.name(s$file), ': no valid corrs'))
  } else {
    ## Some of these corrs may be NA if the firing rate is low, but they
    ## should be safely ignored in the plot.
  
    if (is.null(main)) {
      main = paste(basename(s$file), "dt:", s$corr$dt)
    }
  
    if (is.null(show.pts))
      show.pts <- s$NCells < 100

    if (show.pts) {
      
      plot.default(dists, corrs, xlab=xlabel, ##log=log,
                   ylab=ylabel, bty="n",
                   main=main, col=dot.col,
                   ...)
    } else {
      ## set the ylim to a sensible default.
      upper.pts <- s$corr$corr.id.means[,"mean"] + s$corr$corr.id.means[,"sd"]
      ylim <- c(0.001, max(upper.pts, na.rm=TRUE)) #sd could be NA
      
      plot.default(dists, corrs, xlab=xlabel, type='n',
                   ylab=ylabel, bty="n",
                   main=main, ylim=ylim,
                   ...)
      show.ci <- TRUE                   #better show something.
    }


    if (identify) {
      labels1 <- outer(seq(1, s$NCells), seq(1,s$NCells), FUN="paste")
      labs <- labels1[which(upper.tri(labels1))]
      identify(dists, corrs, labels=labs)
    }

    if (show.ci) 
      .plotCI(s$corr$corr.id.means[,"mid"], s$corr$corr.id.means[,"mean"],
             s$corr$corr.id.means[,"sd"],
             xlab=xlabel, ylab=ylabel,
             pch=19, add=TRUE)
    if (show.fit) 
      .corr.do.fit(s$corr$corr.id,plot=TRUE)

    if (!is.null(s$corr$method) && show.method) {
      title(sub = s$corr$method)
    }
  }
}

.write.corr.indexes <- function(s, file=NULL) {
  ## Write out the correlation index values to a CSV file for
  ## further processing.
  ncells = s$NCells
  op = matrix(0, nrow= (ncells*(ncells-1)/2), ncol=4)

  colnames(op) = c("unit.i", "unit.j", "distance (um)", "corr index")
  d = s$corr$dists                      #distance matrix
  c = s$corr$corr.indexes               #correlation matrix
  n=1;
  for (j in 1:(ncells-1)) {
    for (i in (j+1):ncells) {
      op[n,1] = j;
      op[n,2] = i;
      op[n,3] = d[j,i];
      op[n,4] = c[j,i];
      n=n+1
    }
  }

  if (is.null(file)) {
    file <- paste(get.project.plate.name(s$file),"_corrs.csv",sep="")
    cat(sprintf("Writing correlations to %s\n", file))
  }
  write.csv(op, file=file, row.names=FALSE)
  ## Return the file as well in case we want it.
  invisible(op)
}

.plot.corr.index.fit <- function(s, ...) {
  ### Show the correlation indexes and then the fit.
  .plot.corr.index(s, identify=FALSE,col="red", log="")
  .plotCI(s$corr.id.means[,1], s$corr.id.means[,2], s$corr.id.means[,3],
         xlab="distance", ylab="correlation index", 
         pch=19, add=TRUE)
  .corr.do.fit(s$corr.id,plot=TRUE)
}


## TODO: this loop is quite slow for large numbers of spike trains,
## and probably could be rewritten in C for speed.

.make.corr.indexes <- function(spikes, dt, min.rate=0) {
  ## Return the correlation index values for each pair of spikes.
  ## The matrix returned is upper triangular.
  ## SPIKES should be a list of length N, N is the number of electrodes.
  ## "dt" is the maximum time for seeing whether two spikes are coincident.
  ## This is defined in the 1991 Meister paper.
  ## If MIN.RATE is >0, use the electrode iff the firing rate is above
  ## MIN.RATE.
  
  n <- length(spikes)
  if (n == 1) {
    ## If only one spike train, cannot compute the cross-corr indexes.
    0;
  } else {
    Tmax <- max(unlist(spikes))           #time of last spike.
    Tmin <- min(unlist(spikes))           #time of first spike.

    no.minimum <- isTRUE(all.equal(min.rate, 0))

    if (!no.minimum) {
      ## precompute rates, and find which electrodes are okay.
      rates <- sapply(spikes, length) / (Tmax - Tmin)
      rates.ok <- rates > min.rate
      .printf('Rejecting %d electrodes with firing rate below %.3f Hz\n',
             n-sum(rates.ok), min.rate)
    }
    
    corrs <- array(0, dim=c(n,n))
    for ( a in 1:(n-1)) {
      n1 <- length(spikes[[a]])
      for (b in (a+1):n) {
        n2 <- length(spikes[[b]])
        if ( no.minimum || (rates.ok[a] && rates.ok[b])) {
          val <- as.double(.count.nab(spikes[[a]], spikes[[b]],dt) *
                           (Tmax-Tmin)) /
                             (as.double(n1) * n2 * (2*dt))
          if (is.na(val)) {
            stop(sprintf('.make.corr.indexes: NA generated for pair %d %d',
                         a, b))
          }
        } else {
          ## one of the electrodes was below min firing rate.
          val <- NA
        }
        corrs[a,b] <- val
      }
    }

    corrs
  }
}


.make.corr.indexes2 <- function(spikes, dt, min.rate=0) {
  ## New version using the C routine for corr indexing.
  ## Return the correlation index values for each pair of spikes.
  ## The matrix returned is upper triangular.
  ## SPIKES should be a list of length N, N is the number of electrodes.
  ## "dt" is the maximum time for seeing whether two spikes are coincident.
  ## This is defined in the 1991 Meister paper.
  ## If MIN.RATE is >0, use the electrode iff the firing rate is above
  ## MIN.RATE.
  
  n <- length(spikes)
  if (n == 1) {
    ## If only one spike train, cannot compute the cross-corr indexes.
    0;
  } else {
    Tmax <- max(unlist(spikes))           #time of last spike.
    Tmin <- min(unlist(spikes))           #time of first spike.
    
    no.minimum <- isTRUE(all.equal(min.rate, 0))

    if (!no.minimum) {
      ## precompute rates, and find which electrodes are okay.
      rates <- sapply(spikes, length) / (Tmax - Tmin)
      rates.ok <- rates > min.rate
      .printf('Rejecting %d electrodes with firing rate below %.3f Hz\n',
             n-sum(rates.ok), min.rate)
    } else {
      rates.ok <- rep(0, n)             #need to pass to C anyway...
    }
    
    ## corrs <- array(0, dim=c(n,n))

    ## create one long vector of spikes.
    all.spikes <- unlist(spikes)
    nspikes <- sapply(spikes, length)
    duration <- Tmax - Tmin
    
    first.spike <- c(0, cumsum(nspikes)[-n])
    z <- .C("count_overlap_arr",
            as.double(all.spikes),
            as.integer(n),
            as.integer(nspikes),
            as.integer(first.spike),
            as.integer(rates.ok),
            as.integer(no.minimum),
            as.double(duration),
            as.double(dt),
            res = double(n*n))

    ## return the result.
    array(z$res, dim=c(n,n))
  }
}





## ?? This function not used ?? 2009-10-06
## corr.index.means <- function(x) {
##   ## Compute the mean,sd correlation index at each given distance.
##   dists <- x$dists[which(upper.tri(x$dists))]
##   corrs <- x$corr.indexes[which(upper.tri(x$corr.indexes))]

##   dists.uniq <- unique(dists)
##   num.dists <- length(dists.uniq)       #num of  different distances.

##   ##print(dists.uniq)
##   ## create 4-D array to store results.  Each row stores the
##   ## distance, mean corr, sd, and num of values at that distance.

##   res <- array(0,  dim=c(num.dists,4))
##   colnames(res) <- c("dist","mean corr", "sd", "n")
  
##   i <- 1

##   for (d in dists.uniq) {
##     ## find all correlations for pairs within 0.01um of given distance.
##     cs <- corrs[ which(abs(dists-d)<0.01)]
##     corrs.mean <- mean(cs)
##     corrs.sd   <- sd(cs)
##     res[i,] <- c(d, corrs.mean, corrs.sd, length(cs))
##     i <- 1+i
##   }

##   res
## }

.corr.get.means <- function(id, mid) {
  ## mid contains the mid-point of each bin.
  data.by.bin = split(id[,"corr"], id[,"dist.bin"])
  bins.found = as.integer(names(data.by.bin)) #assume sorted?
  mids = mid[bins.found]
  means = sapply(data.by.bin, mean)
  sds = sapply(data.by.bin, sd)
  n = sapply(data.by.bin, length)
  cbind(mid=mids, mean=means, sd=sds, n=n)
}

.corr.get.means.old2 <- function(id) {
  ## Compute the mean,sd of the correlation index at each distance.
  ## id is the array of [n,2] values.  Each row is [d,i].
  ## where d is the distance and i is the correlation.
  ## Returns a matrix.

  dist = id[,1]
  corr = id[,2]

  ## split does the hard work here, of dividing up the correlation
  ## values into those that share the same values of distance.
  
  l = split(corr, dist) 
  m = sapply(l, mean)
  s = sapply(l, sd)
  n = sapply(l, length)
  res = cbind(dist=as.double(names(l)), mean=m, sd=s, n=n)
  rownames(res) <- NULL
  res
}

.corr.get.means.old <- function(id) {
  ## Compute the mean,sd of the correlation index at each distance.
  ## id is the array of [n,2] values.  Each row is [d,i].
  ## where d is the distance and i is the correlation.
  ## Returns a matrix.
  
  corr.get.means.helper <- function(x) {
    ## Helper function to create  mean and sd of one set of distances.
    ## Need to check that the correlation index is not NA.
    ## X is the distance that we are currently processing.
    indexes <- which(( id[,1] == x) & !is.na(id[,2]))
    c(x, mean(id[indexes,2]), sd(id[indexes,2]), length(indexes))
    ##c(x, median(id[indexes,2]), mad(id[indexes,2]), length(indexes))
  }
  
  d.uniq <- sort(unique(id[,1]))
  means <- t(sapply(d.uniq, corr.get.means.helper))
  colnames(means) <- c("dist", "mean", "sd", "n")
  means
}

.corr.get.means.check <- function() {
  ## Simple function to check our new and old functions are the same.

  ## For each test, generate synthetic data and check for same results.
  for (i  in 1:3) {
    n = 1000^2                           #for n cells, need n^2 entries.
    d = sample(1000, n, replace=T)       #imagine 1000 different distances.
    i = runif(n, max=200)
    id = cbind(d, i)
    print(system.time( m1 <- .corr.get.means.old(id)))
    print(system.time( m2 <- .corr.get.means(id)))
    stopifnot(all.equal(m1,m2))
  }
}

.corr.do.fit <- function(id, plot=TRUE, show.ci=FALSE, ...) {
  ## Do the fit to the exponential and optionally plot it.  Any
  ## correlation index of zero is removed, since we cannot take the
  ## log of zero.  Hopefully there won't be too many of these.
  ## If SHOW.CI is true, do the fit with 95% confidence intervals.

  y.zero <- which(id[,2]==0)
  if (length(y.zero)>0) {
    id <- id[-y.zero,]
    warning(paste("removing", length(y.zero),"zero entries"))
  }
  x <- id[,1]
  y.log <- log(id[,2])
  fit <- lm(y.log ~ x)
  if (show.ci) {
    ## TODO: why is 850 hard-coded in here?
    expt.new <- data.frame(x = seq(0, 850, 10))  #range of x for predictions.
    expt.clim <- predict(fit, expt.new, interval="confidence")
  }
  if (plot)  {
    if (show.ci) {
      ## Confidence intervals will show mean, so don't need
      ## to do both matlines and curve.
      matlines(expt.new$x, exp(expt.clim), lty=c(1,2,2),
               col="black")
    } else {
      curve(exp(fit$coeff[1]) * exp(x * fit$coeff[2]), add = TRUE,
            from=0, ...)
    }
  }
  fit
}

.corr.check.fit <- function() {
  ## Simple test routine to see that the exponential fits are okay.
  a <- 40; b <- .01
  x <- seq(from=1, to=500, by=20)
  y <- a*exp(-b*x) + (2*rnorm(length(x)))
  plot(x,y, log="y")
  fit <- .corr.do.fit( cbind(x,y), col='red')
  
  ## should be similar to (a,b)
  print(exp(fit$coefficients))
}

.my.upper <- function (x,diag=FALSE) {
  ## Return the upper triangular elements of a matrix on a
  ## column-by-column basis.
  ## e.g. .my.upper(matrix(1:9, nrow=3), diag=TRUE).
  ## returns >>1 4 5 7 8 9<<
  if (is.matrix(x)) {
   x[ which(upper.tri(x,diag))]
  } else {
    stop(paste(deparse(substitute(x)),"is not a matrix"))
  }
}

.my.lower <- function (x,diag=FALSE) {
  ## Return the lower triangular elements of a matrix on a
  ## column-by-column basis.
  ## e.g. .my.lower(matrix(1:9, nrow=3), diag=TRUE).
  ## returns >>1 2 3 5 6 9<<
  if (is.matrix(x)) {
   x[ which(lower.tri(x,diag))]
  } else {
    stop(paste(deparse(substitute(x)),"is not a matrix"))
  }
}


