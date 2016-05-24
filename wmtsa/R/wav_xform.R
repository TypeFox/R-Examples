################################################
## WMTSA package wavelet transform functionality
##
##  Functions:
##
##    wavBestBasis
##    wavDWPT
##    wavDWT
##    wavMODWPT
##    wavMODWT
##    wavPacketBasis
##    wavPacketIndices
##
##  Constructor Functions and methods:
##
##    wavCWT
##
##      as.matrix.wavCWT
##      plot.wavCWT
##      print.wavCWT
##
##    wavTransform
##
##      [.wavTransform
##      [<-.wavTransform
##      [[.wavTransform
##      as.matrix.wavTransform
##      boxplot.wavTransform
##      eda.plot.wavTransform
##      plot.wavTransform
##      plot.wavTransform.crystal
##      print.wavTransform
##      print.summary.wavTransform
##      wavStackPlot.wavTransform
##      summary.wavTransform
##      reconstruct.wavTransform
##
################################################

######################################################
##
##    wavCWT
##
##      as.matrix.wavCWT
##      plot.wavCWT
##      print.wavCWT
##
######################################################

###
# wavCWT (constructor)
###

"wavCWT" <- function(x, scale.range=deltat(x) * c(1, length(x)), n.scale=100,
  wavelet="gaussian2", shift=5, variance=1){

  # check input argument types and lengths
  checkVectorType(scale.range,"numeric")
  checkScalarType(n.scale,"integer")
  checkScalarType(wavelet,"character")
  checkScalarType(shift,"numeric")
  checkScalarType(variance,"numeric")
  checkRange(n.scale, c(1,Inf))

  # obtain series name
  series.name <- deparse(substitute(x))

  if (length(scale.range) != 2)
    stop("scale.range must be a two-element numeric vector")
  if (variance <= 0)
    stop("variance input must be positive")

  # obtain sampling interval
  sampling.interval <- deltat(x)

  # create a vector of log2-spaced scales over the specified range of scale
  octave <- logb(scale.range, 2)
  scale  <- ifelse1(n.scale > 1, 2^c(octave[1] + seq(0, n.scale - 2) * diff(octave) /
    (floor(n.scale) - 1), octave[2]), scale.range[1])

  # project the scale vector onto a uniform grid of sampling interval width
  scale <- unique(round(scale / sampling.interval) * sampling.interval)
  n.scale <- length(scale)

  # check scale range
  if (abs(min(scale) - sampling.interval) > .Machine$double.eps)
    stop("Minimum scale must be greater than or equal to sampling interval ",
      "of the time series")

  # map the time vector
  if (inherits(x, "signalSeries"))
  {
    times <- as(x@positions,"numeric")
    x <- x@data
  }
  else
  {
  
    times <- time(x)
    x <- as.vector(x)
  }

  # ensure that x is a double vector
  storage.mode(x) <- "double"

  # map the wavelet filter and corresponding argument
  gauss1 <- c("gaussian1", "gauss1")
  gauss2 <- c("gaussian2", "gauss2", "mexican hat", "sombrero")
  supported.wavelets <- c("haar", gauss1, gauss2, "morlet")
  wavelet <- match.arg(lowerCase(wavelet), supported.wavelets)

  # map the filter type to MUTILS type
  # 4: gaussian1
  # 5: gaussian2, sombrero, mexican hat
  # 6: morlet
  # 7: haar
  filter <- mutilsFilterTypeContinuous(wavelet)

  if (filter == 4)
  {
    filter.arg <- sqrt(variance)
    wavelet    <- "gaussian1"
  }
  else if (filter == 5)
  {
    filter.arg <- sqrt(variance)
    wavelet    <- "gaussian2"
  }
  else if (filter == 6)
  {
    filter.arg <- shift
    wavelet    <- "morlet"
  }
  else if (filter == 7)
  {
    filter.arg <- 0.0
    wavelet    <- "haar"

    # in the case of the Haar, the Euler-Macluarin approximation
    # to the DFT of the wavelet filter is not defined for non-integer
    # multiples of the sampling interval. therefore, we coerce the
    # scales accordingly in this case
    scale <- sampling.interval * unique(round(scale / sampling.interval))
  }
  else
  {
    stop("Unsupported filter type")
  }
  
  # calculate the CWT
  z <- itCall("RS_wavelets_transform_continuous_wavelet",
    as.numeric(x),
    as.numeric(sampling.interval),
    as.integer(filter),
    as.numeric(filter.arg),
    as.numeric(scale))
    #
    #COPY=rep(FALSE, 5),
    #CLASSES=c("numeric", "numeric", "integer", "numeric", "numeric"),
    #PACKAGE="ifultools")

  # if the impulse response of the waveleyt filter is purely real
  # then transform CWT coefficients to purely real as well
  if (wavelet != "morlet")
    z <- Re(z)

  # assign attributes
  attr(z, "scale")       <- scale
  attr(z, "time")        <- as.vector(times)
  attr(z, "wavelet")     <- wavelet
  attr(z, "series")      <- x
  attr(z, "sampling.interval") <- sampling.interval
  attr(z, "series.name") <- series.name
  attr(z, "n.sample")    <- length(x)
  attr(z, "n.scale")     <- n.scale
  attr(z, "filter.arg")  <- filter.arg

  oldClass(z) <- "wavCWT"

  z
}

###
# as.matrix.wavCWT
###

"as.matrix.wavCWT" <- function(x, mode="any", names=TRUE, ...)
{
  checkScalarType(mode,"character")
  checkScalarType(names,"logical")
  y <- x
  attributes(y) <- NULL
  dims <- dim(x)
  nrow <- dims[1]
  ncol <- dims[2]
  dimnms <- dimnames(x)
  y <- matrix(as.vector(y, mode=mode), nrow=nrow, ncol=ncol)
  if(names)
    dimnames(y) <- dimnms
  y
}

###
# plot.wavCWT
###

"plot.wavCWT" <- function(x, xlab=NULL, ylab=NULL, logxy="y",
  power.stretch=0.5, phase=FALSE, series=FALSE, series.ylab="",zoom=NULL,
  type="image", grid.size=100, add=FALSE,  theta=120, phi=30, ...)
{
  # define local functions
  "imageScale" <- function(x, power.stretch=0, stretch.fun=NULL, ...){
    if(is.null(stretch.fun)) {
      if(power.stretch == 0)
	return(log(abs(x) + 1))
	else return((abs(x))^power.stretch)
    }
    else stretch.fun(x)
  }

  # set plot layout
  if (!add){
   frame()
    # store image plot settings for possible
    # overlay of WTMM skeleton
    if (is.element(type, "image")){
      on.exit(par(fig=fig))
      on.exit(par(usr=usr))
    }

    old.plt <- splitplot(1,1,1)
    on.exit(par(old.plt))
  }

  zz <- as.matrix(x)
  if(is.complex(x))
    zz <- ifelse1(phase, Arg(zz), Mod(zz))

  # obtain attributes
  xatt     <- attributes(x)
  times    <- xatt$time
  scales   <- xatt$scale
  x.series <- xatt$series
  n.sample <- xatt$n.sample
  n.scale  <- length(scales)

  x.series.units <- ifelse1(is(x.series,"signalSeries"),
    x.series@units.position, NULL)

  if (is.null(xlab))
    xlab <- ifelse1(is.null(x.series.units), "Time", x.series.units)
  if (is.null(ylab))
    ylab <- ifelse1(is.null(x.series.units), "log2(scale)",
      paste("Scale: log2(", x.series.units, ")"))
  type <- match.arg(type, c("image", "persp"))

  if (type == "image"){

    plt <- par("plt")

    if (series){
      old.fig <- par(fig=c(0,1,0,0.9))
      on.exit(par(old.fig))
    }

    data        <- scaleZoom(times, scales, zz, zoom=zoom, logxy=logxy)
    series.time <- data$x
    itime       <- data$ix

    image(data$x, data$y, imageScale(data$z, power.stretch=power.stretch), ...,
      xlab=xlab, ylab=ylab, axes=TRUE)

    # store plot information for future possible
    # overlay of WTMM skeleton
    usr <- par("usr")
    fig <- par("fig")

    if (series){

      mai <- par("mai")
      old.fig <- par(fig=c(fig[1:2], 0.8, 0.95), new=TRUE)
      on.exit(par(old.fig))

      plot(series.time, x.series[itime], type="l", col="blue", axes=TRUE,
        xlim=range(series.time), ylab=series.ylab,xlab="",xaxt="n", xaxs="i")
      if (prod(sign(range(x.series))) < 0)
	lines(par("usr")[1:2], rep(0,2), lty="dashed")
    }
  }
  else if (type == "persp"){

    iscale <- seq(1, n.scale, length=min(grid.size, n.scale))
    itime  <- seq(1, n.sample, length=min(grid.size, n.sample))

    data <- scaleZoom(times[itime], scales[iscale], zz[itime,iscale], zoom=NULL, logxy=logxy)

    persp(data$x, data$y, data$z,
	      xlab=xlab, ylab=ylab, zlab="CWT Modulus", theta=theta, phi=phi, ...)
   }
  else
    stop("Unsupported plot type")

  invisible(NULL)
}

###
# print.wavCWT
###

"print.wavCWT" <- function(x, justify="left", sep=":", ...)
{
  # obtain attributes
  xatt       <- attributes(x)
  wavelet    <- lowerCase(xatt$wavelet)
  filter.arg <- xatt$filter.arg
  name       <- xatt$series.name
  scale      <- xatt$scale
  n.scale    <- length(scale)
  series     <- xatt$series
  sampling.interval <- xatt$sampling.interval

  # pretty print strings
  waveletstr <- switch(wavelet,
    haar      = "Haar",
    gaussian1 = "Gaussian (first derivative)",
    gaussian2 = "Mexican Hat (Gaussian, second derivative)",
    morlet    = "Morlet",
    wavelet)

  wavelet <- match.arg(wavelet, c("gaussian1", "gaussian2", "morlet"))

  filtval1 <- filtval2 <- NULL

  if (any(wavelet == c("gaussian1", "gaussian2"))){

    filtcat1 <- "Wavelet variance"
    filtval1 <- filter.arg^2

    if (is(series, "signalSeries")){
      units.time <- series@units.position
      if (length(units.time) > 0)
        filtval1 <- paste(filtval1, " (", units.time, ")", sep="")
    }
  }

  if (wavelet == "morlet"){
    filtcat2 <- "Wavelet frequency shift"
    filtval2 <- paste(filter.arg, "(rad/sec)")
  }

  scale.range <- range(scale)

  main <- paste("Continuous Wavelet Transform of", name)

  z <- list(
    "Wavelet"=waveletstr,
    "Wavelet variance"=filtval1,
    "Wavelet frequency shift"=filtval2,
    "Length of series"=length(series),
    "Sampling interval"=sampling.interval,
    "Number of scales"=n.scale,
    "Range of scales"=paste(scale.range[1], "to", scale.range[2])
  )

  prettyPrintList(z, header=main, sep=sep, justify=justify, ...)

  invisible(x)
}

###
# wavBestBasis
###

"wavBestBasis" <- function(costs)
{
  # check input arguments
  checkVectorType(costs,"numeric")

  # define local functions
  flatIndex <- function(j,n, index.base=1) 2^j - 1 + n + index.base

  # costs is a list containing the DWPT costs
  # in W(0,0), W(1,0), W(1,1), W(2,0), ..., W(J,2^J-1) order
  # costs is a list containing the DWPT costs
  n.level <- ilogb(length(costs), base=2)

  # flatten costs into vector
  costs <- unlist(costs)

  if (length(costs) != 2^(n.level + 1) - 1)
    stop("number of costs not representative of a wavelet packet transform")

  # allocate memory for output
  basis <- vector("integer", length(costs))

  # initialize basis vector
  last.level <- flatIndex(n.level,0)
  for (i in seq(length(basis)))
	  basis[i] <- ifelse1(i < last.level, 0, 1)

  # loop through the decomposition levels
  # from the penultimate to zero
  for (j in seq(n.level - 1, 0)){

	  # loop through the oscillation indices (local
	  # node indices) from left to right

	  for (b in seq(0,2^j-1)){

	    # define parent and left child index
	    b1 <- flatIndex(j,b)
	    b2 <- flatIndex(j+1, 2*b)

      # compare cost of parent to the sum of the cost
      # of the children. if the parent
      # has a lower cost than the children, then the
      # parent node is 'marked' as a possible best basis
      # node. otherwise, the cost of the parent is replaced by
      # the sum of the children nodes.
      parent   <- costs[b1]
      children <- costs[b2] + costs[b2+1]

	    if (parent < children){

        # mark the parent node as a possible best basis crsytsal
        basis[b1] <- 1

        # zero out the subtree below current parent node.
	      # this eliminates all of the children as possible
	      # candidates for a best basis
	      for (jj in seq(j+1, n.level)){

	        K <- 2^(jj-j)

	        for(k in seq(b*K, K*(b+1) - 1)){
	          basis[flatIndex(jj,k)] <- 0
	        }
	      }

	    }
	    else{

	      # eliminate the parent as a possible best basis crystal candidate
	      basis[b1] <- 0

	      # set the cost of the parent node to the sum of the children's costs,
	      # a lower cost than that of the parent
        costs[b1] <- children
      }
   }
 }

 x     <- which(basis > 0)-1
 level <- ilogb(x + 1, base=2)
 osc   <- x - 2^level + 1

 # repack costs into list
 costs <- lapply(seq(0,n.level),function(j,costs) {costs[seq(2^j, 2^(j+1)-1)]}, costs=costs)

 return(list(level=level,	osc=osc, costs=costs))
}

###
# wavDWPT
##

"wavDWPT" <- function(x, wavelet="s8", n.levels=ilogb(length(x), base=2),
  position=list(from=1,by=1,units=character()), units=character(),
  title.data=character(), documentation=character())
{
  if (is(x, "wavTransform"))
    return(x)

  checkScalarType(n.levels,"integer")
  checkScalarType(wavelet,"character")

  # obtain the wavelet and scaling filters
  filters <- wavDaubechies(wavelet=wavelet, normalized=FALSE)

  # set the default number of decomposition levels
  # to be the largest scale at which there exists
  # at least one interior (non-boundary) wavelet
  # coefficient
  N <- length(x)
  L <- length(filters$wavelet)

  if (is.null(n.levels))
    n.levels <- max(wavMaxLevel(n.taps=L, n.sample=N, xform="dwpt"), 1)

  if (L > N)
	  stop("Filter length must exceed or be equal to the length of time series")

  # get the series name
  series.name <- deparse(substitute(x))

  # convert data to signalSeries class
  if (!length(title.data))
    title.data <- series.name

  x <- create.signalSeries(x, position=position, units=units,
    title.data=title.data, documentation=documentation)

  # calculate the transform
  data <- itCall("RS_wavelets_transform_packet",
    as.numeric(x),
    list(filters$wavelet,filters$scaling),
    as.integer(n.levels))
    #
    #COPY=rep(FALSE,3),
    #CLASSES=c("numeric","list","numeric"),
    #PACKAGE="ifultools")

  # convert each object from a matrix to a vector
  data <- lapply(data, as.vector)

  # create crystal names
  crystals <- NULL

  for (j in seq(0,n.levels))
    crystals <- c(crystals, paste("w",j,".", seq(0,2^j-1),sep=""))

  names(data) <- crystals

  # name each atom in each crystal
  for (j in seq(data))
    names(data[[j]]) <-
      paste(names(data[j]), "(", seq(0,length(data[[j]])-1), ")",sep="")

  # create dictionary
  dictionary <- wavDictionary(wavelet=wavelet, dual=FALSE, decimate=FALSE,
	  n.sample=length(x),	attr.x=NULL, n.levels=n.levels,	boundary="periodic",
	  conv=TRUE, filters=filters,	fast=TRUE, is.complex=FALSE)

  wavTransform(data, x, as.integer(n.levels), dictionary, FALSE, "dwpt")
}

###
# wavDWT
###

"wavDWT" <- function(x, n.levels=ilogb(length(x), base=2),
  wavelet="s8", position=list(from=1,by=1,units=character()), units=character(),
  title.data=character(), documentation=character(), keep.series=FALSE)
{
  if (is(x, "wavTransform"))
    return(x)

  checkScalarType(n.levels,"integer")
  checkScalarType(wavelet,"character")

  # get the series name from the parent
  # (or from this function if this is called as a top level function)
  series.name <- deparseText(substitute(x))

  # convert data to signalSeries class
  if (!length(title.data))
    title.data <- series.name

  series <- create.signalSeries(ifelse1(keep.series, x, NULL),
    position=position, units=units,
    title.data=title.data,
    documentation=documentation)

  # map filter number and obtain length
  filters <- wavDaubechies(wavelet=wavelet, normalized=FALSE)

  # perform transform
  data <- lapply(
    itCall("RS_wavelets_transform_discrete_wavelet_convolution",
      as.numeric(x),
      list(filters$wavelet,filters$scaling),
      n.levels),
    as.vector)

  # develop names of the crystals
  crystals <- c(paste("d", 1:n.levels, sep=""), paste("s",n.levels,sep=""))

  if (length(data) == (n.levels + 2)){

    crystals <- c(crystals,"extra")

    # create names for extra coefficients
    n     <- length(x)
    i     <- 1
    extra <- rep("",length(data[[n.levels + 2]]))

    for (j in 1:n.levels){

      if (n %% 2){ # is odd
	    extra[i] <- paste("s",j-1,sep="");
	    i <- i + 1;
	    n <- n - 1
      }
      n <- n / 2
    }

    names(data[[n.levels + 2]]) <- extra
  }

  # name the crystals in the return list
  names(data) <- crystals

  # name each atom in each crystal
  for (j in seq(data)){
    names(data[[j]]) <-
      paste(names(data[j]), "(", seq(0,length(data[[j]])-1), ")",sep="")
  }

  # obtain filters to store for future use
  # use un-normalized daubechies filters
  filters <- wavDaubechies(wavelet=wavelet, normalized=FALSE)

  # create dictionary
  dict <- wavDictionary(wavelet=wavelet,
    dual=FALSE, decimate=FALSE, n.sample=length(x),
    attr.x=NULL, n.levels =as.integer(n.levels),
    boundary="periodic", conv=TRUE, filters=filters,
    fast=TRUE, is.complex=FALSE)

  # use wavTransform constructor to create output object
  wavTransform(data=data, series=series, n.levels=as.integer(n.levels),
    dictionary=dict, shifted=FALSE, xform="dwt")
}

###
# wavMODWPT
###

"wavMODWPT" <- function(x, wavelet="s8", n.levels=ilogb(length(x), base=2),
  position=list(from=1,by=1,units=character()), units=character(),
  title.data=character(), documentation=character())
{
  if (is(x, "wavTransform"))
    return(x)

  checkScalarType(n.levels,"integer")
  checkScalarType(wavelet,"character")

  # get the series name from the parent
  # (or from this function if this is called as a top level function)
  series.name <- deparseText(substitute(x))

  # convert data to signalSeries class
  if (!length(title.data))
    title.data <- series.name

  # obtain the wavelet and scaling filters
  filters <- wavDaubechies(wavelet=wavelet, normalized=TRUE)

  # set the default number of decomposition levels
  # to be the largest scale at which there exists
  # at least one interior (non-boundary) wavelet
  # coefficient
  N <- length(x)
  L <- length(filters$wavelet)

  if (is.null(n.levels))
    n.levels <- max(wavMaxLevel(n.taps=L, n.sample=N, xform="modwpt"), 1)

  if (L > N)
    stop("Filter length must exceed or be equal to the length of time series")

  x <- create.signalSeries(x, position=position, units=units,
		title.data=title.data, documentation=documentation)

  n.sample <- length(x)

  # calculate the MODWPT
  data <- itCall("RS_wavelets_transform_maximum_overlap_packet",
    as(x,"numeric"),
    list(filters$wavelet, filters$scaling),
    as.integer(n.levels))
    #COPY=rep(FALSE,3),
    #CLASSES=c("numeric","list","numeric"),
    #PACKAGE="ifultools")

  # convert each object from a matrix to a vector
  data <- lapply(data, as.vector)

  # create crystal names
  crystals <- NULL

  for (j in seq(0,n.levels))
   crystals <- c(crystals, paste("w",j,".", seq(0,2^j-1),sep=""))

  names(data) <- crystals

  # name each atom in each crystal
  for (j in seq(data))
    names(data[[j]]) <- paste(names(data[j]), "(", seq(0,n.sample-1), ")",sep="")

  # create dictionary
  dictionary <- wavDictionary(wavelet=wavelet,
    dual=FALSE, decimate=FALSE, n.sample=length(x),
    attr.x=NULL, n.levels =as.integer(n.levels),
    boundary="periodic", conv=TRUE, filters=filters,
    fast=TRUE, is.complex=FALSE)

  wavTransform(data, x, as.integer(n.levels), dictionary, FALSE, "modwpt")
}

###
# wavMODWT
###

"wavMODWT" <- function(x, wavelet="s8", n.levels=ilogb(length(x), base=2),
	 position=list(from=1,by=1,units=character()), units=character(),
	 title.data=character(), documentation=character(), keep.series=FALSE)
{
  if (is(x, "wavTransform"))
    return(x)

  checkScalarType(n.levels,"integer")
  checkScalarType(wavelet,"character")

  # define local functions
  "sortCrystals" <- function(x, reverse=FALSE)
  {
    # check for proper class
    if (!(class(x) == "character"))
      stop("Input must be of class character")

    # develop local sort function
    localsort <- function(x) x[order(as.numeric(substring(x,2)))]

    # divide the string into wavelet and scaling crystals
    wave.crystals <- x[grep("d",x)]
    scal.crystals <- x[grep("s",x)]

    # sort each individually and recombine
    sorted.crystals <- c(localsort(wave.crystals), localsort(scal.crystals))

    # reverse if requested
    if (!reverse)
      return(sorted.crystals)
    else
      return(rev(sorted.crystals))
  }

  if (is(x, "wavTransform")){
     if (x$xform == "modwt")
       return(x)
     else
       stop ("Cannot convert input transform to MODWT")
  }

  # obtain the wavelet and scaling filters
  filters <- wavDaubechies(wavelet=wavelet, normalized=TRUE)

  # set the default number of decomposition levels
  # to be the largest scale at which there exists
  # at least one interior (non-boundary) wavelet
  # coefficient
  if (is.null(n.levels)){
    L <- length(filters$wavelet)
    N <- length(x)

    if (N > L)
      n.levels <- ilogb(((N - 1) / (L - 1) + 1), base=2)
    else if (N > 0)
      n.levels <- 1
    else
      stop("length of time series must be positive")
  }

  # get the series name from the parent
  # (or from this function if this is called as a top level function)
  series.name <- deparseText(substitute(x))

  # convert data to signalSeries class
  if (!length(title.data))
    title.data <- series.name

  series <- create.signalSeries(ifelse1(keep.series, x, NULL),
    position=position, units=units, title.data = title.data,
    documentation=documentation)

  # call the modwt
  data <- lapply(
    itCall("RS_wavelets_transform_maximum_overlap",
      as.numeric(x),
      list(filters$wavelet,filters$scaling),
      n.levels),
      #COPY=rep(FALSE,3),
      #CLASSES=c("numeric","list","numeric"),
      #PACKAGE="ifultools"),
    as.vector)

  # create dictionary
  dict <- wavDictionary(wavelet=wavelet,
    dual=FALSE, decimate=FALSE, n.sample=length(x),
    attr.x=NULL, n.levels =as.integer(n.levels),
    boundary="periodic", conv=TRUE, filters=filters,
    fast=TRUE, is.complex=FALSE)

  # create crystal names
  if((J <- n.levels) > 0)
    crystals <- c(paste("d", 1:J, sep=""), paste("s",J, sep=""))
  else
    crystals <- "s0"

  # split rows of wavelet coefficients matrix into a list
  names(data) <- crystals

  # name each atom in each crystal
  for (j in seq(data)){
    names(data[[j]]) <- paste(names(data[j]), "(",
      seq(0,length(data[[j]])-1), ")",sep="")
  }

  # use wavTransform constructor to create output object
  wavTransform(data=data, series=series, n.levels=as.integer(n.levels),
    dictionary=dict, shifted=FALSE, xform="modwt")
}

##############################################################################
##
##  Constructor wavTransform:
##
##      [.wavTransform
##      [<-.wavTransform
##      [[.wavTransform
##      as.matrix.wavTransform
##      boxplot.wavTransform
##      eda.plot.wavTransform
##      plot.wavTransform
##      plot.wavTransform.crystal
##      print.wavTransform
##      print.summary.wavTransform
##      wavStackPlot.wavTransform
##      summary.wavTransform
##      reconstruct.wavTransform
##
##############################################################################

###
# Constructor: wavTransform
###

"wavTransform" <- function(data, series, n.levels, dictionary, shifted, xform)
{
  # check input arguments
  if (!is.list(data))
    stop("data must be a list")
  if (!all(unlist(lapply(data, isVectorAtomic))))
    stop("data list must contain vector objects as defined by isVectorAtomic()")
  if (!isVectorAtomic(series) && !is(series, "signalSeries"))
    stop("series must be a vector as defined by isVectorAtomic() or ",
      "an object of class \"signalSeries\"")
  if (!is.integer(n.levels) || length(n.levels) > 1 || n.levels < 1)
    stop("n.levels must be a positive integer scalar")
  if (!is(dictionary, "wavDictionary"))
    stop("dictionary must be an object of class \"wavDictionary\"")
  if (!is.logical(shifted) || length(shifted) > 1)
    stop("shifted be a logical scalar")
  if (!is.character(xform) || length(xform) > 1)
    stop("xform be a character string scalar")

  # construct object list
  z <- list(
    data       = data,
    scales     = 2^(0:(n.levels-1)),
    series     = series,
    dictionary = dictionary,
    shifted    = shifted,
    xform      = xform)

  oldClass(z) <- "wavTransform"

  return(z)
}

###
# [.wavTransform
###

"[.wavTransform" <- function(x, ...)
{
  # define local functions
  "levelNodes" <- function(x, level){

  	n.level <- x$dictionary$n.levels

  	if (level < 1 || level > n.level)
  	  stop(paste("Level must be on the interal [1,",n.level,"]",sep=""))

  	if (is.element(x$xform, c("dwt","modwt"))){
  	  index <- ifelse1(level == n.level,
  	    paste(c("d","s"), level, sep=""), level)
  	}
  	else if (is.element(x$xform, c("dwpt","modpwt"))){
  	  osc   <- seq(0,(2^level)-1)
	  index <- (2^level) + osc
  	}
  	else
  	  stop("Data extraction method not supported for transform type ", x$xform)

  	index
  }

  # check input arguments
  index <- ..1
  if (!is.character(index) && is.numeric(index))
    index <- as.integer(index)
  if (!is.character(index) && !is.integer(index))
    stop("Index must be an object of class character or integer")

  # check to see if a level has been specified
  index <- ifelse1(is.element("level", names(list(...))),
    levelNodes(x, list(...)$level), index)

  # retain only the specified indices
  x$data <- x$data[index]

  x
}

###
# [<-.wavTransform
###

"[<-.wavTransform" <- function(x, i, ..., value){

  by.crystal <- FALSE

  # obtain number of coefficients in each crystal
  ncoeff <- sapply(x$data,length)

  # obtain number of replacement coeffs
  nrep <- length(value)

  # check to see if replacement object is a list.
  if (is.list(value))
    by.crystal <- TRUE

  if (!by.crystal){

    # if length of replacement value vector
    # is less than the number of coefficients
    # in each scale, then replicate the last
    # value in the replacement vector accordingly

    # make the substitutions
    for (crystal in i){

      crys <- names(x$data[[crystal]])

      replacement <- value

      # if the replacement is too long for current crystal, then truncate it.
      # if it is too short, replicate last coefficient to make up the remainder
      if (nrep > ncoeff[crystal])
        replacement <- replacement[1:ncoeff[crystal]]
      else if (nrep < ncoeff[crystal])
        replacement <- c(replacement, rep(replacement[nrep], ncoeff[crystal] - nrep))

      replacement       <- unlist(replacement)
      x$data[[crystal]] <- replacement

      names(x$data[[crystal]]) <- crys
    }
  }
  else{

    # check to make sure that the number of
    # replacement values equals the number of
    # crystals to be replaced
    ilen <- length(i)

    if (nrep > ilen)
      value <- value[1:ilen]
    else if (nrep < ilen)
      value <- c(value,rep(value[length(value)],ilen-nrep))

    # replace the crystals with the respective values.
    # we use the match() function to find the correct
    # value to replace with in the case where the
    # crystals are defined as strings like c("d2","d1")
    # for example
    for (crystal in i){

      replacementData <- value[[match(crystal,i)]]

      len <- length(replacementData)

      if (len == ncoeff[crystal])
        x$data[[crystal]] <- replacementData
      else
        x$data[[crystal]] <- c(replacementData, rep(replacementData[len],
          ncoeff[crystal] - len))
    }
  }

  x
}

###
# [[.wavTransform
###

"[[.wavTransform" <- function(x, ...)
{
	index <- ..1
  if (!is.character(index) && is.numeric(index))
    index <- as.integer(index)

  if (!is.character(index) && !is.integer(index))
    stop("Index must be an object of class character or integer")

  index <- index[1]

  if (is.integer(index)){

  	n.crystal <- length(x$data)

  	if (index < 1 || index > n.crystal)
  	  stop("Index must be on the interval [1,", n.crystal, "]", sep="")
  }
  else if (!is.element(index, names(x$data)))
    stop("Crystal name does not exist")

  x$data[[index]]
}

###
# as.matrix.wavTransform
###

"as.matrix.wavTransform" <- function(x, names=TRUE, ...)
{
  y <- x$data
  attributes(y) <- NULL
  if (!names)
    return(unlist(x$data, use.names=FALSE))

  nms  <- unlist(lapply(x$data,names),use.names=FALSE)
  data <- unlist(x$data,use.names=FALSE)
  names(data) <- nms
  as.matrix(data)
}

###
# boxplot.wavTransform
###

"boxplot.wavTransform" <- function (x, xlab="Crystal",
  ylab="Distribution Statistics", ...)
{
  crystals <- rev(names(x$data))
  boxplot(x$data[crystals],labels=crystals, xlab=xlab, ylab=ylab, ...)
  invisible(NULL)
}

###
# eda.plot.wavTransform
###

"eda.plot.wavTransform" <- function(x, data=TRUE, n.top=NULL, cex.main=1, mex=.5,
  x.legend=NULL, y.legend=NULL, gap=.13, ...)
{
  old.plt <- splitplot(2,2,1, gap=gap)
  on.exit(par(old.plt))

  dict <- x$dictionary
  J <- dict$n.levels
  if(J==0)
    cat("Trivial Transform\n")
  else{
    plot(wavShift(x), add=TRUE, cex.main=cex.main, adj=1)
    splitplot(2,2,2, gap=gap)
    plot(x, type="energy", plot.pie=FALSE, add=TRUE)
    splitplot(2,2,3, gap=gap)
    boxplot(x, cex=1, ...)
    title("Crystal Distribution", cex.main=cex.main, adj=1, line=0.5)
    splitplot(2,2,4, gap=gap)
    plot(x, type="energy", plot.bar=FALSE, add=TRUE, cex.main=cex.main)
  }
  invisible(NULL)
}

###
# plot.wavTransform
###

"plot.wavTransform" <- function(x, type="h", plot.bar=TRUE, plot.pie=TRUE,
  add=FALSE, vgap=0.05, grid=TRUE, grid.lty="dashed", border=TRUE, cex.main=1, ...)
{
  "energy.plot" <- function(x, tit="Energy Distribution", plot.bar=TRUE, plot.pie=TRUE, add=FALSE, cex.main=0.7)
  {
     if (!plot.pie && !plot.bar)
      plot.bar <- TRUE

     same.plot <- as.logical(plot.bar && plot.pie && !add)

     if (!add){
       old.plt <- ifelse1(same.plot, splitplot(2,1,1), splitplot(1,1,1))
       on.exit(old.plt)
     }

     energy <- unlist(lapply(x$data,function(x) sum(x^2)))
     colors <- seq(along=energy)+1
     nms    <- names(x$data)

     if (plot.bar){
       barplot(energy, names=nms, col=colors)
       title(tit, cex.main=cex.main, adj=1)
     }

     if (plot.pie){

       if (same.plot)
         splitplot(2,1,2)

       # for the pie chart, combine energies that are
       # small so as not to muddle the pie chart
       tol  <- 0.1
       bad  <- which(as.logical(energy / max(energy) < tol))
       good <- which(as.logical(energy / max(energy) >= tol))

       if (length(bad) > 0){
         energy  <- c(energy[good], sum(energy[bad]))
         pie.str <- c(nms[good], "Other")
         colors  <- c(colors[good], 200)
       }
       else{
         energy  <- energy[good]
         pie.str <- nms[good]
         colors  <- colors[good]
       }

       pie.str <- paste(paste(pie.str, ": ",sep=""),
          round(energy/sum(energy)*100,2),"%",sep="")

       pie(energy, labels=pie.str, radius=0.95,
          clockwise=FALSE, col = colors)
 
       if (!same.plot)
         title(tit, cex.main=cex.main, adj=1)
     }
  }

  "packet.plot" <- function(x, vgap=0.05, grid=TRUE, grid.lty="dashed",
    border=TRUE, add=FALSE, cex.main=0.7, ...){

    singleplot <- as.logical(!add)

    if (singleplot){
      frame()
      old.plt <- splitplot(1,1,1)
      par(c(old.plt, list(new=FALSE)))
    }

    # turn off clipping so that
    # labels outside of plot appear
    old.xpd <- par(xpd=NA)
    on.exit(par(old.xpd))

    if (is.complex(x))
      stop("cannot plot complex data")

    dict <- x$dictionary
    n.sample <- dict$n.sample
    n.levels <- dict$n.levels

    # if the number of crystals is less than
    # that in a full transform, the user has
    # extracted a subset and thus only
    # issue a wavStackPlot
    if (length(x$data) < (2^(n.levels+1) - 1)){
      wavStackPlot.default(x$data)
      return(NULL)
    }

    tt <- seq(0, n.sample - 1)

    # obtain all crystal names
    crystals <- names(x$data)

    # obtain min and max of data
    yrange <- range(unlist(x$data))
    ygap <- abs(diff(yrange)) * vgap
    ymin <- yrange[1] - ygap
    ymax <- yrange[2] + ygap
    dy   <- abs(ymax - ymin)
    dy2  <- dy/2

    # plot the data
    plot(tt, seq(0,(n.levels+1)*dy,length=length(tt)), type="n", axes=FALSE,
	    ylim=c(0, (n.levels+1)*dy), xlab="", ylab="", ...)

    start <- 0
    em    <- par("cxy")[1L]
    left  <- par("usr")[1L] - em

    for (j in seq(0,n.levels)){

	    nodes <- seq(2^j)
      ybias <- (n.levels - j) * dy
	    ydata <- unlist(x$data[nodes+start]) - ymin
	    xdata <- seq(0,n.sample-1,length=length(ydata))
	    col   <- ifelse1(!j, "blue", "black")

	    lines(xdata, ydata + ybias, col=col, ...)
      text(left, ybias + dy2, paste("Level",j))

      # add y-axis limits to first plot
	    if (!j){
	      right <- par("usr")[2]
	      text(right, (n.levels+1)*dy, paste(round(ymax,2)))
	      text(right, (n.levels)*dy, paste(round(ymin,2)))
	    }

	    start <- start + length(nodes)
    }

    J <- n.levels
    tt.delta <- tt[2]-tt[1]
    tt.start <- tt[1]
    tt.end <- tt[n.sample]

    if (border){
      segments(tt.start, c(0, J+1)*dy, tt.end, c(0, J+1)*dy, lty=1, lwd=3)
      segments(c(tt.start, tt.end), 0, c(tt.start, tt.end), (J+1)*dy, lty=1, lwd=3)
    }

    if (grid){

      # draw horizontal grid lines
      segments(tt.start, (0:(J+1))*dy, tt.end, (0:(J+1)) * dy, lty=grid.lty)
      segments(c(tt.start, tt.end), 0, c(tt.start, tt.end), (J+1)*dy, lty=grid.lty)

      # draw vertical grid lines
      for(i in 1:J){
	    tti <- tt.start+(tt.end-tt.start)*seq(1, 2^i-1, 2)/2^i
	    segments(tti, (J-i+1)*dy, tti, 0, lty=grid.lty)
      }
    }
  }

  if (is.element(type, "energy")){
    energy.plot(x, plot.bar=plot.bar, plot.pie=plot.pie, add=add, cex.main=cex.main, ...)
    invisible(NULL)
  }
  else{

      if (is.element(x$xform, c("dwt","modwt")))
      invisible(wavStackPlot(x, cex.main=cex.main, same.scale=TRUE, type=type, add=add, ...))
  	else if (is.element(x$xform, c("dwpt","modwpt")))
  	  invisible(packet.plot(x, vgap=0.05, grid=TRUE,
  	  grid.lty=grid.lty, border=TRUE, add=add, cex.main=cex.main, ...))
  }
}

###
# plot.wavTransform.crystal
###

"plot.wavTransform.crystal"  <- function(x, ...)
  invisible(wavStackPlot(x, ...))

###
# print.wavTransform
###

"print.wavTransform" <- function(x, justify="left", sep=":", ...)
{
  series.name <- x$series@title
  dict        <- x$dictionary
  crystals    <- names(x$data)

  xform <- switch(x$xform,
	  dwt    = "Discrete Wavelet Transform",
	  modwt  = "Maximal Overlap Discrete Wavelet Transform",
	  dwpt   = "Discrete Wavelet Packet Transform",
	  modwpt = "Maximal Overlap Discrete Wavelet Packet Transform",
	  character(0))

  many.crystals <- as.logical(length(crystals) > 12)

  main <- paste(xform, "of", series.name)
  z <- c(as.list(dict),
    list("Zero phase shifted"=x$shifted,
      "Crystals"=ifelse1(many.crystals,"", crystals)))
  prettyPrintList(z, header=main, sep=sep, justify=justify, ...)

  if(many.crystals){
	  crys.mat <- matrix(crystals[1:12],nrow=4,ncol=3,byrow=TRUE)
	  dimnames(crys.mat) <- list(rep("",nrow(crys.mat)),rep("",ncol(crys.mat)))
	  print(crys.mat, quote=FALSE)
	  cat(paste(" ... (",length(crystals)," bases)\n",sep=""))
  }

  invisible(x)
}

###
# print.summary.wavTransform
###

"print.summary.wavTransform" <- function(x, digits=max(2, .Options$digits - 4), ...)
{
  cat("\n")
  print(round(oldUnclass(x$smat), digits=digits), ...)
  cat("\nEnergy Distribution:\n")
  print(round(x$Edist[1:3,], digits=digits), ...)
  cat("\n")
  invisible(x)
}

###
# wavStackPlot.wavTransform
###

"wavStackPlot.wavTransform" <- function(x, data=TRUE, same.scale=TRUE,
  title.=NULL, xlab=NULL, cex.main=0.7, zeroline=TRUE, times=NULL, add=FALSE, adj=1, ...)
{
  # save plot parameters and restore upon exit
  if (!add){
   frame()
   old.plt <- splitplot(1,1,1)
   on.exit(par(old.plt))
  }

  series         <- x$series
  series.omitted <- as.logical(length(series@data) == 0)
  position.units <- series@units.position
  data.units     <- series@units

  if (is.null(xlab)){
    if (is.null(series@title)) xlab <- "Position"
    else{
      if (length(position.units) > 0) xlab <- position.units
      else xlab <- "Position"
    }
  }

  crystals    <- names(x$data)
  series.name <- wavTitle(x)
  wavelet     <- x$dictionary$wavelet

  # remove extra scaling coefficients
  iextra <- which(crystals == "extra")
  if (length(iextra > 0))
    crystals <- crystals[-iextra]

  # create title
  if (is.null(title.))
	title. <- paste(upperCase(x$xform), "of",
		series.name, "using", wavelet, "filters")

  # combine original sequence with transform coefficients
  if (series.omitted)
    y <- x$data[crystals]
  else{

    y <- c(list(series@data), x$data[crystals])

    if (nchar(series.name) > 8)
	  series.name <- paste(substring(series.name,1,5), "...",sep="")

    if (length(data.units) > 0)
	  series.name <- paste(series.name,"\n(",data.units,")",sep="")

    names(y[1]) <- series.name
  }

  if (x$shifted){

    zerophaseshifts <- switch(x$xform,
	    dwt=wavIndex(x)$shift$dwt[crystals],
	    modwt=wavIndex(x)$shift$modwt[crystals],
	    NA)

    if (series.omitted)
      names(y) <- paste(crystals,"{",zerophaseshifts,"}",sep="")
    else
      names(y) <- c(names(y[1]), paste(crystals,"{",zerophaseshifts,"}",sep=""))
  }

  # obtain series information
  if (is.null(times)){

    if (series.omitted) times <- seq(length=x$dictionary$n.sample)
    else times <- as(series@positions,"numeric")
  }

  wavStackPlot.default(y, same.scale=same.scale, y.axis=TRUE,
	 cex.main=cex.main, zeroline=zeroline, times=times,
	 bars=FALSE, ...)

  title(main=title., cex.main=cex.main, adj=adj)
  mtext(text=xlab, side=1, line=4, cex=cex.main)
  invisible(NULL)
}

###
# summary.wavTransform
###

"summary.wavTransform" <- function(object, ...)
{

  # calculate statistics on each crystal of the data
  crystals   <- names(object$data)

  measures <- c("Min", "1Q", "Median",
		"3Q", "Max", "Mean",
		"SD", "Var","MAD", "Energy %")

  smat <- matrix(0, length(crystals), length(measures))

  E <- unlist(lapply(crystals,
		     function(crystal, data) sum(data[[crystal]]^2),
		     data=object$data))

  smat.list <- lapply(crystals,
    function(crystal,E,x){
      xi     <- oldUnclass(x[[crystal]])
      quantx <- quantile(xi, c(0, .25, .5, .75, 1))
      meanx  <- mean(xi)
      varx   <- var(xi)
      stdx   <- sqrt(varx)
      madx   <- mad(xi)
      enormx <- sum(xi^2)/E*100
        c(quantx,meanx,stdx,varx,madx,enormx)
    },
    E =sum(E),
    x =object$data)

  smat <- do.call("rbind",smat.list)

  dimnames(smat) <- list(crystals, measures)

  # calculate the energy distribution
  energy.pcent <- c(1:5, 10, 15, 20, 25)

  Edist <- matrix(0, 3, length(energy.pcent)+1)

  dimnames(Edist) <- list(c("Energy %", "|coeffs|", "#coeffs"),
    c("1st", paste(energy.pcent, "%", sep="")))

  allE <- oldUnclass(unlist(object$data[crystals]))^2
  E <- sum(allE)		# total energy
  ii <- c(1, ceiling(length(allE)*energy.pcent/100))
  allE <- rev(sort(allE))
  Edist[1,] <- cumsum(allE/E*100)[ii]
  Edist[2,] <- sqrt(allE[ii])
  Edist[3,] <- ii

  obj  <- list(smat=smat, Edist=Edist)

  oldClass(obj) <- "summary.wavTransform"
  obj
}

###
# reconstruct.wavTransform
###

"reconstruct.wavTransform" <- function(x, indices=0, ...)
{
  dict <- x$dictionary

  if(dict$n.levels == 0)
    return(as.vector(x))

  # obtain the wavelet and scaling filters

  wavelet.filter <- dict$analysis.filter$high
  scaling.filter <- dict$analysis.filter$low

  # call the inverse transform function
  if (x$xform == "dwt"){

    synthesis <- itCall("RS_wavelets_transform_discrete_wavelet_convolution_inverse",
      lapply(x$data,as.vector),
      list(wavelet.filter,scaling.filter))
      #COPY=rep(FALSE, 2),
      #CLASSES=rep("list", 2),
      #PACKAGE="ifultools")
  }
  else if (x$xform == "modwt"){

    synthesis <- as.vector(itCall("RS_wavelets_transform_maximum_overlap_inverse",
      lapply(x$data, function(x) matrix(x, nrow=1)),
      list(wavelet.filter,scaling.filter)))
      #COPY=c(FALSE,FALSE),
      #CLASSES=c("list","list"),
      #PACKAGE="ifultools"))
  }
  else if (x$xform == "dwpt"){

    indices <- wavPacketIndices(indices)$flat
    basis   <- wavPacketBasis(x, indices=indices)

    if (is.null(basis$extra)){

      # create faux extra atoms and levelmap vectors
      data     <- basis
      nextra   <- as.integer(0)
      atoms    <- matrix(as.integer(0))
      levelmap <- matrix(as.integer(0))
    }
    else{
      data     <- basis$data
      atoms    <- matrix(as.integer(basis$extra$atoms), nrow=1)
      levelmap <- matrix(as.integer(basis$extra$levelmap), nrow=1)
      nextra   <- as.integer(length(atoms))
    }

    data <- lapply(data, function(x) matrix(as.double(x), nrow=1))

    synthesis <- itCall("RS_wavelets_transform_packet_inverse",
      data, as.integer(nextra), atoms, levelmap, as.integer(indices),
      list(wavelet.filter,scaling.filter))
      #COPY=rep(FALSE,6),
      #CLASSES=c("list", "integer", rep("matrix",3), "list"),
      #PACKAGE="ifultools")

    synthesis <- as.vector(synthesis)
  }
  else
    stop("Transform is not supported")

  as.vector(synthesis)
}

###
# wavPacketBasis
##

"wavPacketBasis" <- function(x, indices=0)
{
	# Returns the DWPT crystals (in a list) corresponding to the
	# basis specified by the indices vector. The indices
	# are mapped as follows:
	#
	# 0   - original series
	# 1:2 - W(1,0), W(1,1)    (i.e., all level 1 crystals)
	# 3:6 - W(2,0),...,W(2,3) (i.e., all level 2 crystals)
	#
	# and so on. If the indices do not form a basis, an error is issued
	if (!is(x,"wavTransform") && !is.element(x$xform, "dwpt"))
    stop("Input must be an object of class wavTransform and",
      " must be a DWPT")

  zz <- lapply(x$data, as.matrix)

  z <- itCall("RS_wavelets_transform_packet_basis",
    zz, matrix(as.integer(indices)))
    #COPY=rep(FALSE,2), #CLASSES=c("list", "matrix"),
    #PACKAGE="ifultools")

  z <- lapply(z, as.vector)

  indices <- wavPacketIndices(indices)$flat + 1

  n.crystal <- length(indices)
  any.extra <- length (z) == n.crystal + 2
  nms       <- names(x$data)[ indices ]

  data <- z[ seq(n.crystal) ]
  names(data) <- nms

  if (any.extra){

    extra <- z[ seq(n.crystal + 1, n.crystal + 2) ]
    names(extra) <- c("atoms", "levelmap")
    levels <- seq(0, length(extra$levelmap) - 1)

    # name the extra levels vector
    names(extra$levelmap) <- paste("j=", levels, sep="")

    # name the extra atoms
    nms <- NULL
    Nj  <- length(x$data[[1]])

    for (j in levels){

      if (extra$levelmap[ j + 1 ] >= 0){

        osc <- seq(from=0, length=2^j)
        nms <- c(nms, paste("w(", j, ",", osc, ",", Nj, ")", sep=""))
        Nj  <- Nj - 1
      }

      Nj <- Nj / 2
    }

    names(extra$atoms) <- nms

    return(list(data=data, extra=extra))
  }
  else
    return(data)
}

###
# wavPacketIndices
##

"wavPacketIndices" <- function(x, check.basis=TRUE)
{
  checkVectorType(x,"integer")
  if (any(x < 0))
    stop("Indices must be non-negative")

  # define local functions
  flatIndex <- function(j,n, index.base=1) 2^j - 1 + n + index.base

  # obtain level and osc
  level <- ilogb(x+1, base=2)
  osc   <- x - 2^level + 1

  # check for wavelet packet basis
  if (check.basis){

    n.level <- max(level)
    n.crystal <- 2^(n.level + 1) - 1

    basis <- rep(0, n.crystal)

    # wipe out flattened index vector up to lowest input index
    if (min(x) > 0)
      basis[seq(min(x))] <- 1

    for (i in seq(along=x)){

      b1 <- flatIndex(level[i], osc[i])

      # mark the parent node and eliminate all children
      if (basis[b1] == 1)
       stop("Basis violation: there is an overlap in frequency content corresponding ",
         "to the specified wavelet packet indices")

      basis[b1] <- 1
      j <- level[i]
      n <- osc[i]

      # zero out the subtree below current parent node.
	    # this eliminates all of the children as possible
	    # candidates for a best basis
	    if (j < n.level){
	      for (jj in seq(j+1, n.level)){

	        K <- 2^(jj-j)

	        for(k in seq(n*K, K*(n+1) - 1)){

	          b1 <- flatIndex(jj,k)
	          if (basis[b1] == 1)
              stop("Basis violation: there is an overlap in frequency content corresponding ",
                "to the specified wavelet packet indices")
	          basis[b1] <- 1
	        }
	      }
	    }
    }

    last.level <- basis[seq(n.crystal-2^n.level+1, n.crystal)]

    if (!all(last.level)){
    	print(basis)
      stop("Flattened indices do not form a wavelet packet basis: ",
        "normalized frequencies on [0,1/2] not sspanned")
    }
  }

  list(flat=x, level=level, osc=osc)
}
