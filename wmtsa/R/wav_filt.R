################################################
## WMTSA filter functionality
##
##  Functions:
##
##    wavCWTFilters
##    wavShift
##    wavZeroPhase
##
##  Constructor Functions and methods:
##
##    wavDaubechies
##
##      plot.wavDaubechies
##      print.wavDaubechies
##
##    wavGain
##
##      plot.wavGain
##      print.wavGain
##
################################################

###
# wavCWTFilters
###

"wavCWTFilters" <- function(wavelet="Gaussian2", frequency=seq(0, 2*pi, length=1000),
  shift=3, variance=1, times=NULL)
{
  checkScalarType(wavelet,"character")
  checkVectorType(frequency,"numeric")
  checkScalarType(shift,"numeric")
  checkScalarType(variance,"numeric")

  wavelet <- match.arg(lowerCase(wavelet),
    c("haar","gaussian1", "gaussian2","sombrero","mexican hat","morlet"))

  # define local functions
  "haar" <- function(x){

    y    <- rep(0, length(x))
    ineg <- which(x > -1 & x <= 0)
    ipos <- which(x > 0 & x <= 1)

    y[ineg] <- - 1 / sqrt(2)
    y[ipos] <- 1 / sqrt(2)

    y
  }

  "gaussian1" <- function(x, s=1)
    (sqrt(2) * x * exp(-x^2 / (2 * s^2)) / s^(3 / 2) / pi^(1/ 4))

  "gaussian2" <- function(x, s=1)
    (2 * (1 - x^2 / s^2) * exp(-x^2 / (2 * s^2)) / sqrt(3 * s) / pi^(1/ 4))

  "morlet" <- function(x, w0=5){
    i <- complex(imaginary=1, real=0)
    C <- pi^(-1 / 4) / sqrt(1 - 4/sqrt(3)*exp(-w0^2 / 4) + sqrt(2) * exp(-w0^2 / 2))
    (C * exp(-i * w0 * x) * (exp(-x^2 / 2) - sqrt(2) * exp(- w0^2 / 4) * exp(-x^2)))
  }

  # map the filter type to MUTILS type
  # 4: gaussian1
  # 5: gaussian2, sombrero, mexican hat
  # 6: morlet
  # 7: haar
  filter  <- mutilsFilterTypeContinuous(wavelet)

  if (!is.null(times)){

    return(switch(filter-3,
      gaussian1(times, sqrt(variance)),
      gaussian2(times, sqrt(variance)),
      morlet(times, shift),
      haar(times)))
  }

  funcarg <- switch(filter-3, sqrt(variance), sqrt(variance), shift, 0.0)

  z <- itCall("RS_wavelets_filters_continuous",
    as.integer(filter),
    as.numeric(funcarg),
    as.numeric(frequency))
    #
    #COPY=rep(FALSE,3),
    #CLASSES=c("integer",rep("numeric",2)),
    #PACKAGE="ifultools")

  z <- as.vector(z)
  attr(z, "frequency") <- frequency
  attr(z, "wavelet")   <- wavelet

  z
}

##############################################################################
##
##  Constructor wavDaubechies:
##
##    plot.wavDaubechies
##    print.wavDaubechies
##
##############################################################################

###
# Constructor: wavDaubechies
###

"wavDaubechies" <- function(wavelet="s8", normalized=TRUE)
{
  # map filter number and obtain length
  filter <- mutilsFilterType(wavelet=wavelet)

  # obtain filters
  data <- itCall("RS_wavelets_filters_daubechies",
    as.integer(filter$length), as.integer(filter$type), as.logical(normalized))
    #
    #COPY=rep(FALSE,3),
    #CLASSES=c("integer","integer","logical"),
    #PACKAGE="ifultools")

  family <- switch(filter$type + 1,
   "Extremal Phase","Least Asymmetric","Best Localized","Coiflet")

  filters <- list(
    wavelet    = as.vector(data[[1]]),
    scaling    = as.vector(data[[2]]),
    family     = family,
    length     = filter$length,
    normalized = normalized,
    name       = wavelet)

  oldClass(filters) <- "wavDaubechies"
  filters
}

###
# plot.wavDaubechies
###

"plot.wavDaubechies" <- function(x, type="time", ...)
{
  "taps" <- function(filter, ...){

    L <- length(filter)
    n <- 0:(L-1)
    plot(n, filter, pch=18, ...)
    segments(n, rep(0,L), n, filter)
    abline(h=0, lty=2)
  }

  # save plot parameters and restore upon exit
  old.plt <- splitplot(2,1,1)
  on.exit(par(old.plt))

  if (charmatch(type, "time", 0)){

    taps(x$wavelet, ylab="Wavelet", xlab="")
    title(paste("IMPULSE RESPONSE: ", x$family, ", length=", x$length, sep=""), cex=0.7)
    splitplot(2,1,2)
    taps(x$scaling, ylab="Scaling", xlab="Tap")
  }
  else if (charmatch(type, "gain", 0) || charmatch(type, "phase", 0)){

     Npad <- 1024
     Nhalf <- Npad / 2 - 1
     pad <- rep(0, max(Npad - x$length, x$length))

     H <- fft(c(x$wavelet,pad))[1:Nhalf]
     G <- fft(c(x$scaling,pad))[1:Nhalf]

     if (charmatch(type, "phase", 0)){
        H <- Arg(H) * 180 / pi
        G <- Arg(G) * 180 / pi
        str <- "PHASE"
        ystr <- "(deg)"
     }
     else{
        H <- abs(H)
        G <- abs(G)
        ystr <- "(deg)"
        str <- "GAIN"
     }

    f <- seq(0, 0.5, length=Nhalf)

    plot(f, H, ylab=paste("Wavelet",ystr), xlab="", type="l")
    title(paste("FREQUENCY RESPONSE ", str, ": ", x$family, ", length=", x$length, sep =""), cex=0.7)
    splitplot(2,1,2)
    plot(f, G, xlab="Normalized Frequency", ylab=paste("Scaling",ystr), type="l")
  }
  else
    stop("Plot type not supported")

  invisible(NULL)
}

###
# print.wavDaubechies
###

"print.wavDaubechies" <- function(x, verbose=TRUE, ...)
{
  cat("Filter type:", x$family, "\n");
  cat("Filter length:", x$length, "\n");

  if (x$normalized)
    norm <- "TRUE"
  else
    norm <- "FALSE"

  cat("Filter normalized: ", norm, "\n")

  if (verbose){
    cat("\nWavelet filter:\n")
    print(x$wavelet)
    cat("\n\nScaling filter:\n")
    print(x$scaling)
  }
  invisible(x)
}

######################################################
## S+Fractal wavelet filter sqaured gain functions
##
## Class wavGain
## Constructor function: wavGain
## Methods:
##
##   plot.wavGain
##   print.wavGain
##
######################################################

###
# wavGain
##

"wavGain" <- function(wavelet="s8", n.levels=5, n.fft=1024, normalize=TRUE)
{
  # exclude biorthogonal wavelets: B-spline and V-spline
  if (any(substring(wavelet,1,1) == c("b","v"))){
    warning("Biorthogonal wavelets are currently unsupported in this routine")
    return(NA)
  }

  # map filter number and obtain length
  filter <- wavDaubechies(wavelet=wavelet, normalized=normalize)
  L      <- length(filter$wavelet)
  npad   <- max(L, n.fft)

  # form impulse response of wavelet and scaling filter
  g <- h <- rep(0, npad)
  g[seq(L)] <- filter$scaling
  h[seq(L)] <- filter$wavelet

  # form gain functions for wavelet and scaling filter
  # set as first row of the output matrix
  G <- H <- matrix(0+0i, nrow=n.levels, ncol=npad)
  G[1,] <- fft(g, inverse=FALSE)
  H[1,] <- fft(h, inverse=FALSE)

  # define Fourier frequencies
  df <- 1 / npad
  f  <- seq(0, npad-1) * df

  # form gain functions for remaining levels
  for (j in seq(2,n.levels)){

     fj <- as.integer(((2^(j-1) * f) %% 1) / df) + 1
     G[j,] <- G[1,fj] * G[j-1,]
     H[j,] <- H[1,fj] * G[j-1,]
  }

  # name the matrix components
  levnames    <- list(paste("Level", seq(n.levels)), f)
  dimnames(G) <- levnames
  dimnames(H) <- levnames

  # create return object (list)
  z <- list(gain=list(low=G, high=H),
    sqrgain= list(low=abs(G)^2, high=abs(H)^2),
    frequency=f,
    wavelet=wavelet,
    n.levels=n.levels,
    normalize=normalize)

  # define class
  oldClass(z) <- "wavGain"

  # return the object
  z
}

###
# plot.wavGain
##

"plot.wavGain" <- function(x, cex.main=1, adj=1, ...)
{

  "plot.wavelet.response" <- function(x, frequency, high=TRUE, ylab="", wavestr=NULL, norm=TRUE, cex.main=0.85, adj=1, ...){
	  n.levels <- nrow(x$low)
	  levels <- seq(n.levels)
	  f.band <- 2^(-(levels+1))
	  f.band <- c(f.band, 0.5, 1 - f.band)
	  midoctave <- 3 / 4 * 2^(-levels)
	  midoctave[1] <- 0.53


	  if (high){
	    tit  <- paste("Wavelet Filter",wavestr)
	    data <- x$high
	    data.max <- 1.1 * apply(abs(data),MARGIN=1,max)

	    xtext <- sort(c(midoctave, 1 - rev(midoctave[-1])))
	    ytext <- c(rev(data.max), data.max[-1])
	    text.str <- c(paste(n.levels:1), paste(2:n.levels))

	  }
	  else{
	    tit  <- paste("Scaling Filter",wavestr)
	    data <- x$low
	    data.max <- 1.1 * apply(abs(data),MARGIN=1,max)

	    xtext <- f.band[ c(levels, (n.levels+2):(2*n.levels+1))]
	    ytext <- data.max / 2 * rep(1, length(xtext))
	    text.str <- paste(c(rev(levels),levels))
	  }

	  if (is.complex(data))
	    data <- Mod(data)

	  if (norm)
	    tit <- paste("Normalized",tit)

	  # plot the zeroth level response
	  ylim <- c(-0.1, 1.1)*max(unlist(data))

	  plot(frequency, data[1,], xlab="Normalized Frequency (f)", ylab=ylab,
		  type="l", ylim=ylim, col=1, ...)

	  title(main=tit, cex.main=cex.main, adj=adj, line=0.5)

	  # add text to help identify curves when printed
	  # in black and white

	  if (high){
	     text(xtext, ytext, text.str)
	  }
	  else{
	     text(xtext[levels], ytext[levels], text.str[levels], adj=0)
	     ix <- seq(n.levels+1, length(xtext))
	     text(xtext[ix], ytext[ix], text.str[ix], adj=1)
	  }

	  if (n.levels > 1){

	    if (n.levels == 2){
		 points(frequency, data[2,], type="l")
	    }
	    else{
		for (j in 2:n.levels){
		  points(frequency, data[j,], type="l", col=j)
	 	}
	    }
	  }

	  # add vertical band lines

	  abline(v=f.band, lty=2)

	  invisible(NULL)
	}

  dots <- list(...)
  if (hasArg("sqrgain"))
    sqrgain <- dots$sqrgain
  else
    sqrgain <- TRUE

  norm    <- x$normalize
  wavestr <- paste("(",x$wavelet,")")

  if (!is(x,"wavGain"))
    stop("Input object has incorrect class")
  else{

    if (sqrgain) {
      old.plt <- splitplot(2,1,1)
      ylab <- "Squared Gain"
      data <- x$sqrgain
      plot.wavelet.response(data, x$frequency, high=TRUE,
	      ylab=ylab, wavestr=wavestr, norm=norm, ...)
	    splitplot(2,1,2)
      plot.wavelet.response(data, x$frequency, high=FALSE,
	      ylab=ylab, wavestr=wavestr, norm=norm, ...)
    }
    else{
      old.plt <- splitplot(2,2,1)
      ylab <- "Gain"
      data <- x$gain

      plot.wavelet.response(data, x$frequency, high=TRUE,
	      ylab=ylab, wavestr=wavestr, norm=norm, ...)
        plot.wavelet.response(data, x$frequency, high=FALSE,
	      ylab=ylab, wavestr=wavestr, norm=norm, ...)

	    splitplot(2,2,2)
      plot(x$frequency,Arg(x$gain$high[1,])*180/pi,type="l",
	      xlab ="Frequency",ylab="Phase (deg)")
        title("Level 1",adj=1,cex=0.7)

      splitplot(2,2,3)
      plot(x$frequency,Arg(x$gain$low[1,])*180/pi,type="l",
	      xlab ="Frequency",ylab="Phase (deg)")
      title("Level 1",adj=1,cex=0.7)
    }
  }

  par(c(old.plt, list(new=FALSE)))

  invisible(NULL)
}

###
# print.wavGain
##

"print.wavGain" <- function(x, justify="left", sep=":", ...)
{
  main <- paste("Gain functions for", x$wavelet,"filters")

  z <- list(
    "Number of levels"=x$n.levels,
    "Number of Fourier frequencies"=length(x$frequency),
    "Filters normalized"=x$normalize)

  prettyPrintList(z, header=main, sep=sep, justify=justify, ...)

  invisible(x)
}

###
# wavShift
###

"wavShift" <- function(x)
{
  # obtain class of inout and store for possible later conversion
  wave.class <- class(x)
  original.wave.class <- wave.class

  # obtain shift statistics and crystal data
  if (is(x, "wavTransform")){

    # obtain and check shift state
    shift.state <- x$shifted

    # obtain dictionary, crystal names, and crystal data
    dict     <- x$dictionary
    crystals <- names(x$data)
    data     <- x$data
  }
  else if (is(x,"wavBoundary")){

    shift.state <- x$shifted

    # obtain dictionary, crystal names, and crystal data
    dict     <- x$dictionary
    crystals <- x$crystals
    data     <- x$all[crystals]
  }
  else
    stop("Input class currently unsupported")

  # remove (mo)dwpt crystal "w0.0"
  bad <- which(as.logical(match(crystals,"w0.0")))
  if (length(bad) > 0)
    crystals <- crystals[-bad]

  # produce vector of acceptable wavelet identifier strings
#  supported.wavelets <- c("s","c")

  # obtain wavelet used in the decomoposition
  wavelet <- dict$wavelet

  # check on the type of wavelet
	#  if (all(is.na(charmatch(supported.wavelets, wavelet)))){
	#    stop(cat("Wavelet type is not appropriate for zero phase ",
	#	"shifting operations. Try using symmlet or coiflet ",
	#	"wavelets as an alternative"))
	#  }

  # obtain zero phase shift factors
  zero <- wavZeroPhase(wavelet=dict$wavelet,
	levels=seq(length=dict$n.levels))

  iextra <- which(crystals == "extra")

  if (length(iextra > 0))
    crystals <- crystals[-iextra]

  if (any(wave.class == c("wavTransform","wavBoundary"))){

    shifts <-  switch(x$xform,
	 dwt=zero$dwt[crystals],
	 modwt=zero$modwt[crystals],
	 dwpt=zero$dwpt[crystals],
	 modwpt=zero$modwpt[crystals],
	 NA)
  }

  if (anyMissing(shifts))
     stop("Zero phase shifts could not be calculated")

  # shift the crystals
  obj <- lapply(seq(along=shifts),
    function(j, data, shifts, forward){

      # obtain shift index for current crystal
      v <- shifts[j]

      # negate shift if forward is FALSE
      if (!forward)
        v <- -v

      # permute the coefficient in the current crystal
      rotateVector(data[[j]], v)

    },
    data=data[crystals],
    shifts=shifts,
    forward=!shift.state)

  names(obj) <- crystals

  # into now repack the results into the original object class
  y <- x

  if (original.wave.class == "wavTransform"){

    # add the appropriate boundary and interior labels to
    # the coefficients if supplied and the shift is negated

    y$data[crystals]  <- obj[crystals]
    y$shifted         <- !shift.state
  }
  else if (original.wave.class == "wavBoundary"){

    y$all[crystals]   <- obj[crystals]

    # now replace the the boundary and
    # interior objects as well

    for (j in 1:length(crystals)){

      crystal <- crystals[j]

      data <- y$all[[crystal]]

      boundary.indices <- which(names(data) == "b")
      interior.indices <- which(names(data) == "nb")

      y$boundary[[crystal]]        <- data[boundary.indices]
      names(y$boundary[[crystal]]) <- paste(boundary.indices)
      y$interior[[crystal]]        <- data[interior.indices]
      names(y$interior[[crystal]]) <- paste(interior.indices)

    }

    y$shifted <- !shift.state
  }

  y
}


###
# wavZeroPhase
###

"wavZeroPhase" <- function(wavelet="s8", levels=1:3)
{
  # initialize variables
  L <- ifelse1(lowerCase(wavelet) == "haar", 2,
    as.integer(substring(wavelet,2)))
  level <- unique(sort(levels))

  # create crystal names
  dwt.crystals <- c(paste("d",levels,sep=""),paste("s",levels,sep=""))

  # produce vector of acceptable wavelet identifier strings
  wavelet <- lowerCase(wavelet)
  if (is.element(wavelet, c("d2","s2","haar")))
    wavelet <- "haar"
  if (is.element(wavelet, c("s4","d4")))
    wavelet <- "d4"
  supported.wavelets <- c("haar", paste("d", c(4,6), sep=""),
    paste("s", seq(8,20,by=2), sep=""), paste("c", seq(6,30,by=6), sep=""))

  # check on the type of wavelet.
  # if acceptable, then calculate the amount of
  # circular shift needed to bring the crystals
  # to approximate zero phase. this is relevant only for Daubechies
  # least asymmetric and coiflet filters.
  map <- mutilsFilterType(wavelet)

  if (!is.element(wavelet, supported.wavelets))
    return(list(dwt=NA, modwt=NA, dwpt=NA, modwpt=NA))

  shifts <- itCall("RS_wavelets_filter_zero_phase",
    map$type,
    map$length,
    max(levels))
    #
    #COPY=rep(FALSE,3),
    #CLASSES=rep("integer",3),
    #PACKAGE="ifultools")

  # parse results
  dwtcols     <- c(levels, levels + max(levels))
  shift.dwt   <- shifts[[1]][dwtcols]
  shift.modwt <- shifts[[2]][dwtcols]
  dwpt        <- shifts[[3]]
  modwpt      <- shifts[[4]]

  shift.dwpt   <- NULL
  shift.modwpt <- NULL

  # grab the appropriate shifts and create crystal
  # names for wavelet packet transform crystals
  for (j in levels){

    first <- 2^j - 1
    last  <- first + 2^j - 1

    dwpt.vals   <- dwpt[first:last]
    modwpt.vals <- modwpt[first:last]
    crystals    <- paste("w",j,".",seq(0,last-first),sep="")

    names(dwpt.vals)   <- crystals
    names(modwpt.vals) <- crystals

    shift.dwpt   <- c(shift.dwpt, dwpt.vals)
    shift.modwpt <- c(shift.modwpt, modwpt.vals)
  }

  # name the components in shift
  names(shift.dwt)   <- dwt.crystals
  names(shift.modwt) <- dwt.crystals

  # return the shifts for the DWT , MODWT,
  # DWPT, and MODWPT for all levels
  list(dwt=shift.dwt, modwt=shift.modwt,
	  dwpt=shift.dwpt, modwpt=shift.modwpt)
}
