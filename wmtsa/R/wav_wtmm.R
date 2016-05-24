######################################################
## WMTSA package Holder spectrum functionality
##
##   holderSpectrum
##   wavCWTPeaks
##
## Class: wavCWTTree
## Constructor function: wavCWTTree
## Methods:
##
##   [.wavCWTTree
##   plot.wavCWTTree
##   print.summary.wavCWTTree
##   print.wavCWTTree
##   summary.wavCWTTree
##
######################################################

###
# holderSpectrum
###

"holderSpectrum" <- function(x, n.scale.min=3, fit=lmsreg)
{
  if (!is(x,"wavCWTTree"))
    stop("x must an object of class wavCWTTree")

  n.branch <- length(x)
  ibranch  <- vector(mode="integer", length=n.branch)
  branch   <- seq(n.branch)
  holder.exponent <- holder.time <- vector(mode="numeric", length=n.branch)

  count <- 0

  for (i in branch){

    times    <- x[[i]]$time
    logscale <- logb(x[[i]]$scale, base=2)
    logwtmm  <- logb(x[[i]]$extrema, base=2)

    if (length(logscale) > 1){

      # order data from small to large scale
      if (logscale[2] < logscale[1]){
        times    <- rev(times)
        logscale <- rev(logscale)
        logwtmm  <- rev(logwtmm)
      }

      # obtain the appropriate scaling range
      # choose the range of scale which corresponds
      # to the smallest scales
      breaks <- linearSegmentation(logscale, logwtmm)
      div    <- ifelse1(length(breaks), breaks[1], length(times))

      if (div >= n.scale.min){

        count <- count + 1
        holder.time[count] <- times[1]
        datalist <- list(logwtmm=logwtmm, logscale=logscale, div=div)

        model <- ifelse1(div > 5, fit(logwtmm ~ logscale, subset=seq(div), data=datalist),
          lm(logwtmm ~ logscale, subset=seq(div), data=datalist))

        holder.exponent[count] <- coef(model)["logscale"]

        ibranch[count] <- i
      }
    }
  }

  igood <- seq(count)
  holder.exponent <- holder.exponent[igood] - 0.5
  holder.time     <- holder.time[igood]
  branch          <- branch[ibranch[igood]]

  itime <- order(holder.time)
  list(exponent=holder.exponent[itime], time=holder.time[itime], branch=branch[itime])
}

###
# wavCWTTree
###

"wavCWTTree" <- function(x, n.octave.min=1, tolerance=0.0, type="maxima")
{
  # define local functions
  "WTMM" <- function(x, tolerance=NULL, type="maxima"){

    if (!is(x,"wavCWT"))
      stop("Input object must be of class wavCWT")

    # obtain attributes
    x.attr   <- attributes(x)
    times    <- x.attr$time
    scales   <- x.attr$scale
    n.sample <- x.attr$n.sample
    series   <- x.attr$series

    # check tolerances
    if (is.null(tolerance)){

      # use Donoho and Johnstone universal thresholding
      # tol=sqrt(2 * sigma^2 * log(N)) where sigma^2
      # is the variance of the noise, and N is the length
      # of the original time series. Since we do not
      # know sigma^2 a priori, we use the median absolute deviation
      # the level 1 DWT coefficients using the Haar filter as an
      # approxiation. Finally, we divide by the sqrt(scale) to make
      # the threshold scale dependent

      # noise.variance <- (median(abs((diff(series)[seq(2,n.sample-1,by=2)] / sqrt(2))))/ 0.6745)^2
      # vprint(noise.variance)
      # tolerance <- sqrt(2 * noise.variance * log(n.sample)) / sqrt(scales)
      tolerance <- mad(Mod(x[,1])) / scales
    }

    if (length(tolerance) < length(scales))
      tolerance <- tolerance[1] / sqrt(scales)

    wtmmz <- itCall("RS_wavelets_transform_continuous_wavelet_modulus_maxima",
      as.matrix(x)+0i, tolerance, mutilsTransformPeakType(type))
    #
      #CLASSES=c("matrix","numeric","integer"),
      #COPY=rep(FALSE,3),
      #PACKAGE="ifultools")

    z <- matrix(0, nrow=nrow(x), ncol=ncol(x))
    z[matrix(unlist(wtmmz),ncol=2)+1] <- 1

    z
  }

  "wtmmBranches" <- function(wtmm, extrema.mask, times, scales, span.min=5, gap.max=3, skip=NULL, sampling.interval=1)
  {
    # define normalized scales
    scales <- as.integer(scales / sampling.interval)

    n.scale <- ncol(extrema.mask)
    n.sample <- nrow(extrema.mask)

    if (is.null(scales))
      scales <- 1:n.scale

    iwtmm <- which(extrema.mask[, n.scale] > 0)

    # scan from large to small scale
    iscale <- seq(n.scale-1,1,-1)

    tree <- as.list(iwtmm)
    names(tree) <- iwtmm
    peakStatus <- as.list(rep(0, length(iwtmm)))
    names(peakStatus) <- iwtmm
    orphanRidgeList <- NULL
    orphanRidgeName <- NULL

    n.level <- length(iscale)

    for (j in seq(n.level)){

      iscale.j <- iscale[j]
      scale.j <- scales[iscale.j]

      if (length(iwtmm) == 0){
          iwtmm <- which(extrema.mask[, iscale.j] > 0)
          next
      }

      span <- scale.j * 2 + 1

      if (span < span.min)
        span <- span.min

      remove.j <- selPeak.j <- NULL

      # loop through each point in the current scale's WTMM
      for (k in seq(along=iwtmm)){

        # define search range in the time index
        itime <- iwtmm[k]
        itime.start <- itime - span
        if (itime.start < 1)
          itime.start <- 1
        itime.end <- itime + span
        if (itime.end > n.sample)
          itime.end <- n.sample
        itime.candidates <- which(extrema.mask[itime.start:itime.end, iscale.j] > 0) + itime.start - 1

        if (length(itime.candidates) == 0){

            status.k <- peakStatus[[as.character(itime)]]

            if (status.k > gap.max & scale.j >= 2){
              temp            <- tree[[as.character(itime)]]
              orphanRidgeList <- c(orphanRidgeList, list(temp[1:(length(temp) - status.k)]))
              orphanRidgeName <- c(orphanRidgeName, paste(iscale.j + status.k + 1, itime, sep="_"))
              remove.j        <- c(remove.j, as.character(itime))
              next
            }
            else {
              itime.candidates <- itime
              peakStatus[[as.character(itime)]] <- status.k + 1
            }
        }
        else {
          peakStatus[[as.character(itime)]] <- 0
          if (length(itime.candidates) >= 2)
            itime.candidates <- itime.candidates[which.min(abs(itime.candidates - itime))]

        }

        tree[[as.character(itime)]] <- c(tree[[as.character(itime)]], itime.candidates)
        selPeak.j <- c(selPeak.j, itime.candidates)
      }

      if (length(remove.j) > 0){
        bad.tree   <- which(is.element(names(tree), remove.j))
        tree       <- tree[-bad.tree]
        peakStatus <- peakStatus[-bad.tree]
      }

      dupPeak.j <- unique(selPeak.j[duplicated(selPeak.j)])

      if (length(dupPeak.j) > 0){

          bad.tree <- NULL

          for (dupPeak.jk in dupPeak.j){
            selInd          <- which(selPeak.j == dupPeak.jk)
            selLen          <- sapply(tree[selInd], length)
            bad.tree.jk     <- which.max(selLen)
            bad.tree        <- c(bad.tree, selInd[-bad.tree.jk])
            orphanRidgeList <- c(orphanRidgeList, tree[bad.tree.jk])
            orphanRidgeName <- c(orphanRidgeName, paste(iscale.j, selPeak.j[bad.tree.jk], sep="_"))
          }

          selPeak.j  <- selPeak.j[-bad.tree]
          tree       <- tree[-bad.tree]
          peakStatus <- peakStatus[-bad.tree]
      }
      names(tree) <- selPeak.j
      names(peakStatus) <- selPeak.j

      if (scale.j >= 2){
        maxInd.next      <- which(extrema.mask[, iscale.j] > 0)
        unSelPeak.j      <- maxInd.next[!is.element(maxInd.next, selPeak.j)]
        newPeak.j        <- as.list(unSelPeak.j)
        names(newPeak.j) <- unSelPeak.j
        tree             <- c(tree, newPeak.j)
        iwtmm            <- c(selPeak.j, unSelPeak.j)
        newPeakStatus    <- as.list(rep(0, length(newPeak.j)))
        names(newPeakStatus) <- newPeak.j
        peakStatus       <- c(peakStatus, newPeakStatus)
      }
      else {
        iwtmm <- selPeak.j
      }
    }

    names(tree) <- paste(1, names(tree), sep="_")
    names(orphanRidgeList) <- orphanRidgeName
    tree <- c(tree, orphanRidgeList)
    tree <- lapply(tree, rev)
    tree <- tree[unique(names(tree))]

    tree <- lapply(seq(along=tree), function(i, tree, iscale.min, times, scales, wtmm){
      itime <- tree[[i]]
      iscale <- seq(iscale.min[i], length=length(itime))
      list(itime=itime, iscale=iscale, time=times[itime], scale=scales[iscale], extrema=wtmm[cbind(itime,iscale)])
    },
    tree=tree,
    iscale.min=as.integer(gsub("_.*","",names(tree))),
    times=times,
    scales=scales*sampling.interval,
    wtmm=wtmm)

    # remove any redundant branches
    iflat <- lapply(tree, function(x, nr) (x$iscale-1)*nr + x$itime, nr=nrow(wtmm))

    flatset <- iflat[[1]]
    bad <- NULL

    for (i in seq(2,length(iflat))){

       if (any(is.element(iflat[[i]], flatset)))
         bad <- c(bad, i)
       else
         flatset <- c(flatset, iflat[[i]])
    }

    if (length(bad) > 0)
      tree <- tree[-bad]

    tree
  }

  # obtain attributes
  x.attr   <- attributes(x)
  times    <- x.attr$time
  scales   <- x.attr$scale
  n.sample <- x.attr$n.sample
  sampling.interval <- x.attr$sampling.interval
  border.times <- range(times) + sampling.interval * c(1,-1)

  # locate the the extrema in the CWT matrix
  extrema.mask <- WTMM(x, tolerance=tolerance, type=type)

  if (!identical(dim(x),dim(extrema.mask)))
    stop("Input WTMM dimenions do not match those of the input CWT matrix")

  # develop WTMM tree
  z <- wtmmBranches(ifelse1(is.complex(x), Mod(as.matrix(x)), as.matrix(x)), extrema.mask, times, scales, sampling.interval=sampling.interval)

  # define minimum number of points needed per branch
  n.scale  <- length(scales)
  n.octave <- log2(max(scales) / min(scales))
  n.voice  <- (n.scale - 1) / n.octave
  n.scale.min <- as.integer(n.voice * n.octave.min)

  good <- which(unlist(lapply(z,function(x, n.scale.min) length(x[[1]]) > n.scale.min, n.scale.min=n.scale.min)))
  z <- z[good]
  endtime <- unlist(lapply(z,function(x,iscale) x$itime[iscale], iscale=which.min(scales)))
  isort <- order(endtime)
  z <- z[isort]

  names(z) <- seq(z)

  attr(z, "iendtime")     <- endtime[isort]
  attr(z, "endtime")      <- times[endtime[isort]]
  attr(z, "time")         <- times
  attr(z, "scale")        <- scales
  attr(z, "extrema.mask") <- extrema.mask
  attr(z, "noise")        <- x[,1]
  attr(z, "branch.hist")  <- colSums(extrema.mask*abs(x))
  attr(z, "wavelet")      <- attr(x,"wavelet")
  attr(z, "filter.arg")   <- attr(x,"filter.arg")
  attr(z, "series.name")  <- attr(x,"series.name")
  attr(z, "series")       <- attr(x,"series")
  attr(z, "sampling.interval") <- attr(x,"sampling.interval")

  oldClass(z) <- "wavCWTTree"

  z
}

###
# [.wavCWTTree
###

"[.wavCWTTree" <- function(x, i, ..., time=NULL, range=NULL)
{
  ax    <- attributes(x)
  times <- ax$endtime

  if (!missing(range)){

    if (length(range) == 2){

      i <- which(times >= range[1] & times <= range[2])
    }
    else{
      i <- seq(length(x))
    }
  }
  else if(!missing(time)){

    i <- sort(time)

    min.scale <- median(diff(times))

    itime <- NULL

    for (j in i){

      itime <- c(itime, which(times >= j - min.scale & times <= j + min.scale))
    }

    i <- itime
  }

  i <- i[i <= length(x) & i >= 1]

  z <- oldUnclass(x)[i]

  attributes(z) <- c(attributes(z),
    ax[ setdiff(names(ax), c("names", "dim", "dimnames")) ])

  z
}

###
# plot.wavCWTTree
###

"plot.wavCWTTree" <- function(x, add=FALSE, pch="o",  label=TRUE, log.="y",
  xlab=NULL, ylab=NULL, extrema=FALSE, zoom=NULL,
  fit=FALSE, models= c("lm", "lmsreg", "ltsreg"),
  cex=0.8, col.skip=max(1,round(256/length(x))),  ...)
{
  branches <- names(x)

  if (!length(x))
    stop("No branches found in input object")

  if (fit){

    n.branch <- min(length(x), 4)

    nrow <- 1
    ncol <- 1

    if (n.branch > 1)
      nrow <- 2
    if (n.branch > 2)
      ncol <- 2

    frame()
    gap <- 0.11
    old.plt <- splitplot(nrow,ncol,1,gap=gap)
    on.exit(par(old.plt))

    lwd <- c(2, 1, 2, rep(1,10))
    lty <- c(1, 1, 4, seq(5,length=10))

    # initialize arguments
    slope <- vector(mode="numeric", length=length(models))

    for (i in seq(n.branch)){

      if (i > 1)
        splitplot(nrow,ncol,i,gap=gap)

      x.branch <- log(x[[i]]$scale)
      y.branch <- log(x[[i]]$extrema)

      if (x.branch[length(x.branch)] < x.branch[1]){
        x.branch <- rev(x.branch)
        y.branch <- rev(y.branch)
      }

      breaks <- linearSegmentation(x.branch, y.branch)

      if (is.null(breaks))
        breaks <- length(x.branch)

      plot(x.branch, y.branch, xlab="log(scale)", ylab="log(|EXTREMA|)", type="b",
        pch="o", cex=cex)

      abline(v=x.branch[breaks], lty=2, xpd=FALSE)

      cut      <- seq(breaks[1]-1)
      x.branch <- x.branch[cut]
      y.branch <- y.branch[cut]

      # the lmsreg() and ltsreg() models will return an
      # error if less than 6 points are fit, so restrict
      # the models used accodingly
      if (length(x.branch) > 5)
        admissible.models <- models
      else
        admissible.models <- "lm"

      datalist <- list(x.branch=x.branch, y.branch=y.branch)

      for (imodel in seq(admissible.models)){

        eval(parse(text=paste("fit <- coef(", admissible.models[imodel],
          "(y.branch ~ x.branch,data=datalist))", sep="")))
        abline(fit, lty=lty[imodel], lwd=lwd[imodel], xpd=FALSE)
        slope[imodel] <- fit[2]
      }

      imodel <- seq(admissible.models)

      key.cex <- 0.7

      mdlName <- upperCase(admissible.models)
	  legend("bottomright",
             paste(format(mdlName), "\t: slope=", format(round(slope[imodel], 4))),
    	     lty=lty[imodel],
    	     lwd=lwd[imodel],
    	     cex=key.cex)

      mtext(paste("Branch", branches[i]), adj=1, cex=0.85, line=0.5)
    }

    return(invisible(NULL))
  }

  if (!add)
    frame()

  x.attr <- attributes(x)

  if (extrema){

    pch <- 20
    col <- "blue"

    mask <- which(x.attr$extrema.mask == 1, arr.ind=TRUE)

    data <- scaleZoom(x.attr$time[mask[,1]], x.attr$scale[mask[,2]], zoom=zoom, logxy=log., xy.linked=TRUE,
      xlab="Time", ylab="Scale")
    if (!add)
      plot(data$x, data$y, col=col, pch=pch, xlab=data$xlab, ylab=data$ylab, cex=cex, ...)
    else
      points(data$x, data$y, col=col, pch=pch, cex=cex, ...)

    return(invisible(NULL))
  }

  if (!add){

    times  <- range(x.attr$time)
    scales <- range(x.attr$scale)

    if(is.element(log., "y")){
      scales <- logb(scales, base=2)
      if (is.null(ylab))
        ylab <- "log2(scale)"
    }
    else if (is.null(ylab))
      ylab <- "Scale"

    if(is.element(log., "x")){
      times <- logb(times, base=2)
      if (is.null(xlab))
        xlab <- "log2(time)"
    }
    else if (is.null(xlab))
      xlab <- "Time"

    plot(times, scales, type="n", xlab=xlab, ylab=ylab)
  }

  col <- c("black","red","blue","green","deeppink","orange","cyan","magenta","violet","navy","purple",
     "yellowgreen","orangered","goldenrod","cornflowerblue","plum","steelblue","tomato","pink")

  for (i in seq(along=x)){

    times  <- x[[i]]$time
    scales <- x[[i]]$scale

    if(is.element(log., "y"))
      scales <- logb(scales, base=2)
    if(is.element(log., "x"))
      times <- logb(times, base=2)

    icol <- ((i-1) %% length(col)) + 1

    if (label){

      imaxscale <- order(scales)[length(scales)]
      lines(times[-imaxscale], scales[-imaxscale], col=col[icol], pch=pch, cex=0.5, type="b", ...)
      text(times[imaxscale], scales[imaxscale], as.character(branches[i]), cex=1, col=col[icol])
    }
    else
      lines(times, scales, col=col[icol], lty=1, pch=pch, ...)
  }

  invisible(NULL)
}


###
# print.summary.wavCWTTree
###

"print.summary.wavCWTTree" <- function(x, digits=max(2, .Options$digits - 4), ...)
{
  oldOptions <- options(digits = digits)
  on.exit(options(oldOptions))
  NextMethod("print")
  invisible(x)

#  print(round(data.frame(oldUnclass(x)), digits), ...)
#  invisible(x)
}

###
# print.wavCWTTree
###

"print.wavCWTTree" <- function(x, justify="left", sep=":", ...)
{
  # obtain attributes
  xatt       <- attributes(x)
  times      <- xatt$time
  name       <- xatt$series.name
  scale      <- xatt$scale
  n.scale    <- length(scale)
  series     <- xatt$series
  sampling.interval <- xatt$sampling.interval
  n.branch   <- length(x)
  filter.arg <- xatt$filter.arg

  # pretty print strings
  waveletstr <- "Mexican Hat (Gaussian, second derivative)"
  filtcat1 <- "Wavelet variance"
  filtval1 <- filter.arg^2

  if (is(series, "signalSeries")){
    units.time <- series@units.position
    if (length(units.time) > 0)
      filtval1 <- paste(filtval1, " (", units.time, ")", sep="")
  }

  scale.range <- range(scale)

  main <- paste("Continuous Wavelet Transform Tree of", name)

  z <- list(
    "Wavelet"=waveletstr,
    "Wavelet variance"=filtval1,
    "Length of series"=length(series),
    "Sampling interval"=sampling.interval,
    "Number of scales"=n.scale,
    "Range of scales"=paste(scale.range[1], "to", scale.range[2]),
    "Number of branches"=n.branch
  )

  prettyPrintList(z, header=main, sep=sep, justify=justify, ...)

  invisible(x)
}

###
# summary.wavCWTTree
###

"summary.wavCWTTree" <- function(object, ...)
{
  tmpargs <- lapply(object, function(x){

    scale <- x$scale
    ext <- x$extrema

    c(x$time[which.min(scale)], length(x$time), log2(max(scale) / min(scale)), range(ext), mean(ext), sd(ext), var(ext), mad(ext))})
    
  z <- data.frame(do.call("rbind", tmpargs))

  names(z) <- c("End Time", "Length", "Octaves", "Min", "Max", "Mean", "SD", "Var", "MAD")

  oldClass(z) <- c("summary.wavCWTTree","data.frame")

  z
}

###
# wavCWTPeaks
###

"wavCWTPeaks" <- function(x, snr.min=3, scale.range=NULL, length.min=10,
  noise.span=NULL, noise.fun="quantile", noise.min=NULL)
{
  if (!is(x,"wavCWTTree"))
   stop("Input must be an object of class wavCWTTree")

  xatt <- attributes(x)
  endtimes  <- attr(x, "endtime")
  times     <- attr(x, "time")
  scale     <- attr(x, "scale")
  noise     <- attr(x, "noise")
  wavelet   <- attr(x,"wavelet")
  series    <- attr(x,"series")
  branch.hist <- attr(x, "branch.hist")
  sampling.interval <- abs(diff(times[1:2]))

  if (!is.element(wavelet,"gaussian2"))
    stop("Only CWT developed using the Mexican hat (gaussian2) filter are supported")

  if (is.null(noise.min))
    noise.min <- quantile(abs(attr(x,"noise")), prob=0.05)

  # define default scale range
  if (is.null(scale.range))
    scale.range <- scale[range(which(branch.hist > quantile(branch.hist,prob=0.8)))]

  # define default noise span
  if (is.null(noise.span))
    noise.span <- max(0.01 * diff(range(times)), 5*sampling.interval)

  # obtain local noise estimates around the termination time of each branch
  noise.levels <- unlist(lapply(endtimes, function(x, noise.fun, times, times.range, noise, noise.min, noise.span){
    time.start <- x - noise.span
    if (time.start < times.range[1])
      time.start <- times.range[1]
    time.end <- x + noise.span
    if (time.end < times.range[2])
      time.end <- times.range[2]

    ix <- which(times >= time.start & times <= time.end)
    noise.local <- noise.fun(abs(noise[ix]))
    if (noise.local < noise.min)
      noise.local <- noise.min

    noise.local
   },
   noise.fun=switch(noise.fun, quantile=function(x){quantile(x, probs=0.95)}, sd=sd, mad=function(x){mad(x, center=0)}),
   times=times,
   times.range=range(times),
   noise=noise,
   noise.min=noise.min,
   noise.span=noise.span))

  tmpargs <- lapply(x, function(x) unlist(lapply(x, function(x, imax) x[imax], imax=which.max(x$extrema))))
  peaks <- data.frame(do.call("rbind", tmpargs))
  peaks <- cbind(data.frame(branch=row.names(peaks)), peaks, data.frame(iendtime=attr(x,"iendtime")))
  peak.snr <- peaks[["extrema"]] / noise.levels
  peak.scale <- peaks[["scale"]]

  branch.lengths <- unlist(lapply(x, function(x, scale.range)
    length(which(x$scale >= scale.range[1] & x$scale <= scale.range[2])),
    scale.range=scale.range))

  # prune branches
  #
  #  good.snr    : estimate of SNR at peak value is greater than or equal to the specified snr.min
  #  good.scale  : the scale of the peak is larger than the minimum of the specified scale range
  #  good.length : the length of the branch within the scale.range is greater than or equal to the specified minimum length.min
  #  good.end    : the index of the terminating time of the branch is on the interval (W, N-W), where N is the length of the time
  #                series and W is integer equivalent of 1/4 the length of the noise span or 3, whichever is greater.
  good.snr     <- peak.snr >= snr.min
  good.scale   <- peak.scale >= scale.range[1]
  good.length  <- branch.lengths >= length.min
  iendtime.min <- max(as.integer(noise.span / sampling.interval / 4),3)
  iendtime.max <- length(times) - iendtime.min + 1
  good.end     <- peaks[["iendtime"]] > iendtime.min & peaks[["iendtime"]] < iendtime.max

  peaks <- peaks[which(good.snr & good.scale & good.length & good.end),]
  row.names(peaks) <- as.character(seq(nrow(peaks)))

  z <- list(x=times[peaks$iendtime], y=series[peaks$iendtime])
  attr(z,"peaks")       <- peaks
  attr(z,"snr.min")     <- snr.min
  attr(z,"scale.range") <- scale.range
  attr(z,"length.min")  <- length.min
  attr(z,"noise.span")  <- noise.span
  attr(z,"noise.fun")   <- noise.fun
  attr(z,"noise.min")   <- noise.min
  z
}
