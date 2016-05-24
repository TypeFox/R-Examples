#############################################################
## WMTSA fractionally differenced process functionality
##
##  Functions:
##
##    wavFDPBand
##    wavFDPBlock
##    wavFDPSDF
##    wavFDPTime
##
##  Constructor Functions and methods:
##
##    wavFDP
##
##      eda.plot.wavFDP
##      plot.wavFDP
##      print.wavFDP
##      print.summary.wavFDP
##      summary.wavFDP
##
#############################################################

###
# wavFDPBand
##

"wavFDPBand" <- function(delta=1/4, method="bandpass", scaling=TRUE,
  levels=1:5, n.sample=n.sample <- 2^(max(levels)+1))
{
  # force the levels to be monotonic and uniform if scaling is FALSE
  # i.e. ensure that j=1, ..., J and that no levels are skipped
  if (scaling && length(levels) != max(levels))
      stop("If scaling is TRUE, the levels vector must contain ",
        "indices j=1, 2, ..., max(levels)")

  method <- match.arg(method, c("bandpass", "integration"))

  if (method == "integration"){

    # develop band limits of integration

    lower.limit <- 2^(-(levels + 1))
    upper.limit <- 2^(-levels)
    mult.factor <- 2^(levels + 1)

    # set the process variance to unity for now until we establish
    # an MLE for delta
    variance <- 1

    # calculate the average SDF value at the center frequency for
    # each band, i.e. C'_j for j=1:J
    Cprime <- unlist(lapply(levels,
      function(j, delta, lim.low, lim.high, mult.factor){

        y <- integrate(f=wavFDPSDF, lower=lim.low[j], upper=lim.high[j],
          delta=delta, variance=1, response=NULL)

	      mult.factor[j] * y$value
      },
	    delta    = delta,
	    lim.low  = lower.limit,
      lim.high = upper.limit,
      mult.factor = mult.factor)
    )

    "fdpsdf.bandpass.scaling" <- function(n.sample, delta, variance, Cprime)
    {
        J <- length(Cprime)
        levels <- seq_len(J)
        n.sample * ( variance * gamma(1 - 2 * delta) / gamma(1 - delta)^2 - sum(Cprime / 2^levels) )
    }
        
    # now calculate the average SDF value
    # for the scaling coefficient(s) band
    if (scaling){
      CprimeJP1 <- fdpsdf.bandpass.scaling(n.sample, delta, variance, Cprime)
      Cprime    <- c(Cprime, CprimeJP1 )
    }
  }
  else if (method == "bandpass"){

    if (!scaling)
      n.sample <- -1

    Cprime <- as.vector(itCall("RS_wavelets_fdp_bandpass_variance",
      as.integer(levels),
	    as.numeric(delta),
	    as.integer(n.sample))
    #,
	    #COPY=rep(FALSE,3),
	    #CLASSES=c("integer","numeric","integer"),
            #PACKAGE="ifultools")
        )
  }
  else{
    stop("Unsupported method specified for mid-octave SDF value approximation")
  }

  # return the vector of average SDF bands
  crystals <- paste("d", levels, sep="")
  if (scaling)
    crystals <- c(crystals, paste("s",max(levels),sep=""))

  names(Cprime) <- crystals

  Cprime
}

###
# wavFDPBlock
##

"wavFDPBlock" <- function(x, wavelet="s8", levels=NULL, sdf=NULL,
  boundary=NULL, edof.mode=1,
  estimator="wlse", delta.range=c(-10.0,10.0),
  position=list(from=1,by=1,units=character()), units=character(),
  title.data=character(), documentation=character(), keep.series=FALSE )
{
  # get the series name from the parent
  # (or from this function if this is called as a top level function)
  series.title <- deparse(substitute(x))

  # map estimator
  estimator <- match.arg(estimator, c("wlse", "mle"))

  # map filter number and obtain length
  filter <- mutilsFilterType(wavelet=wavelet)

  # obtain filters (needed for storage in dictionary)
  filters <- wavDaubechies(wavelet=wavelet, normalized=TRUE)

  # store series
  series <- ifelse1(keep.series,
    create.signalSeries(x, position=position, units=units,
      title.data=title.data, documentation=documentation),
    create.signalSeries(NULL, position=position, units=units,
      title.data=series.title, documentation=documentation))

  # map estimator
  estimator <- lowerCase(estimator)
  full.decomposition <- FALSE

  # set default for boundary variable and make
  # sure that it coordinates with the estimator
  if (is.null(boundary))
    boundary <- ifelse1(estimator == "mle", "nonstationary", "unbiased")

  if (estimator == "wlse"){
      estim  <- 1
      imatch <- charmatch(boundary, c("biased", "unbiased"))

      if (is.na(imatch))
	    stop("\nThe weighted least squares estimator requires that\n",
			"the boundary mode be either \"biased\" or \"unbiased\"")

      if (imatch == 1)
        bound <- list(mode=TRUE, description="biased")
      else if (imatch == 2)
        bound <- list(mode=FALSE, description="unbiased")
      else
        stop("Boundary method is unsupported")
  }
  else if (estimator == "mle"){

      estim  <- 0
      imatch <- charmatch(boundary, c("stationary", "nonstationary"))

      if (is.na(imatch))
        stop("\nThe maximum likelihood estimator requires that\n",
          "the boundary mode be either \"stationary\" or \"nonstationary\"")

      if (imatch == 1){
        bound <- list(mode=TRUE, description="stationary FD process model")
        full.decomposition <- TRUE
      }
      else if (imatch == 2)
        bound <- list(mode=FALSE,
          description="stationary-nonstationary FD process model")
      else
        stop("Boundary method is unsupported")

  }
  else
    stop("Unsupported estimator method")

  # if not supplied, define the LEVELS vector
  max.level <- wavMaxLevel(n.taps=filter$length, n.sample=length(x),
	  xform=switch(estimator, mle="modwt", wlse="dwt"))

  if (is.null(levels))
    levels <- ifelse1(length(x) > filter$length,
      ifelse1(max.level > 2, seq(from=2, length=max.level - 1), seq(1, max.level)), 1)

  if (full.decomposition)
    levels <- seq(1, max.level)

  if (max(levels) > max.level)
   stop(paste("\nNumber of levels exceeds maximum:", max.level,
     "\nEither reduce the filter length or maximum decomposition level."))

  names(levels) <- paste("d",levels,sep="")

  # calculate the block-averaged FD model
  # parameter estimates
  if (is.null(sdf) || !length(sdf)){

	  sdf <- (1:2)
	  sdf.empty <- TRUE
  }
  else
    sdf.empty <- FALSE

  data <- itCall("RS_wavelets_fdp_estimator_block",
    as.double(x),
    as.integer(levels),
    as.integer(filter$type),
    as.integer(filter$length),
    as.integer(estim),
    as.logical(bound$mode),
    as.integer(edof.mode),
    as.double(sdf),
    as.double(delta.range))
#,
    #COPY=rep(FALSE,9),
    #CLASSES=c("numeric",rep("integer",4),"logical","integer","numeric","numeric"),
    #PACKAGE="ifultools")

  # remove values outside of delta range. these values
  # are assigned by the C code as a flag for the case
  # where there are not enough interior coefficients
  # in which to form a maximum likelihood estimation
  na.indices <- which(data[[1]] < min(delta.range) | data[[1]] > max(delta.range))
  data <- lapply(data, function(x, na.indices){
    x[na.indices] <- NA
    x}, na.indices=na.indices)

  # create a dictionary
  dictionary <- wavDictionary(wavelet=wavelet,
    dual=FALSE, decimate=FALSE, n.sample=length(x),
    attr.x=NULL, n.levels=length(levels),
    boundary="periodic", conv=TRUE, filters=filters,
    fast=TRUE, is.complex=FALSE)

  # pack results into a list and return
  wavFDP(estimator=estimator,
    delta=as.vector(data[[1]]),
    variance.delta=ifelse1(estimator == "mle", NA, as.vector(data[[2]])),
    innovations.variance=ifelse1(estimator == "wlse", NA, as.vector(data[[3]])),
    delta.range=delta.range,
    dictionary=dictionary,
    levels=as.integer(levels),
    edof.mode=as.integer(edof.mode),
    boundary=bound,
    series=series,
    sdf.method="Integration lookup table",
    type="block")
}

###
# wavFDPSDF
##

"wavFDPSDF" <- function(f, delta=0.45, variance=1, response=NULL)
{
  mult.fact <- ifelse1(!is.null(response),
    approx(response$frequency, response$sqrgain, f, rule=2)$y, 1)

  mult.fact * variance / (2 * sin(pi * f))^(2 * delta)
}

###
# wavFDPTime
##

"wavFDPTime" <- function(x, wavelet="s8", levels=NULL,
  biased=FALSE, estimator="mle",
  dof.order=0, delta.range=c(-10.0,10.0),
  position=list(from=1,by=1,units=character()), units=character(),
  title.data=character(), documentation=character(), keep.series=FALSE)
{
  # get the series name from the parent
  # (or from this function if this is called as a top level function)
  series.title <- deparse(substitute(x))

  # map filter number and obtain length
  filter <- mutilsFilterType(wavelet=wavelet)

  # obtain filters (needed for storage in dictionary)
  filters <- wavDaubechies(wavelet=wavelet, normalized=TRUE)

  # set the default number of decomposition levels
  # to be the largest scale at which there exists
  # at least one interior (non-boundary) wavelet
  # coefficient
  if (is.null(levels)){

    L <- filter$length
    N <- length(x)

    if (!biased){
      if (N > L)
        n.levels <- ilogb(((N - 1) / (L - 1) + 1), base=2)
      else if (N > 0)
        n.levels <- 1
      else stop("length of time series must be positive")
    }
    else
      n.levels <- 2

    levels <- seq(n.levels)
  }

  names(levels) <- paste("d",levels,sep="")

  # store series
  series <- ifelse1(keep.series,
    create.signalSeries(x, position=position, units=units,
      title.data=title.data, documentation=documentation),
    create.signalSeries(NULL, position=position, units=units,
	  title.data=series.title, documentation=documentation))

  # test estimator
  estimator <- match.arg(lowerCase(estimator), c("mle","lse"))
  method    <- as.integer(estimator == "lse")

  # obtain MLEs of FDP paramter delta
  # and corresponding innovation variance

  data <- itCall("RS_wavelets_fdp_estimator_instantaneous",
    as.numeric(x),
    as.integer(levels),
    as.integer(filter$type),
    as.integer(filter$length),
    as.integer(method),
    as.logical(biased),
    as.integer(dof.order),
    as.numeric(delta.range))
    #,
    #COPY=rep(FALSE,8),
    #CLASSES=c("numeric",rep("integer",4),"logical","integer","numeric"),
    #PACKAGE="ifultools")

  # remove values outside of delta range. these values
  # are assigned by the C code as a flag for the case
  # where there are not enough interior coefficients
  # in which to form a maximum likelihood estimation
  na.indices <- which(data[[1]] < min(delta.range) | data[[1]] > max(delta.range))
  data <- lapply(data, function(x, na.indices){
    x[na.indices] <- NA
    x}, na.indices=na.indices)

  # create a dictionary
  dictionary <- wavDictionary(wavelet=wavelet,
    dual=FALSE, decimate=FALSE, n.sample=length(x),
    attr.x=NULL, n.levels=length(levels),
    boundary="periodic", conv=TRUE, filters=filters,
    fast=TRUE, is.complex=FALSE)

  # pack results into a list and return
  boundary <- ifelse1(biased, list(mode=TRUE, description="biased"),
    list(mode=FALSE, description="unbiased"))

  # pack results into a list and return
  wavFDP(estimator=estimator,
    delta=as.vector(data[[1]]),
    variance.delta=as.vector(data[[2]]),
    innovations.variance=as.vector(data[[3]]),
    delta.range=delta.range,
    dictionary=dictionary,
    levels=as.integer(levels),
    edof.mode=as.integer(dof.order),
    boundary=boundary,
    series=series,
    sdf.method="Integration lookup table",
    type="instantaneous")
}

################################################
##
##  Constructor Functions and methods:
##
##    wavFDP
##
##      eda.plot.wavFDP
##      plot.wavFDP
##      print.wavFDP
##      print.summary.wavFDP
##      summary.wavFDP
##
################################################


###
# wavFDP (constructor)
##

"wavFDP" <- function(estimator, delta, variance.delta,
  innovations.variance, delta.range, dictionary, levels,
  edof.mode, boundary, series, sdf.method, type)
{
  # check input arguments
  checkScalarType(estimator, "character")
  checkVectorType(delta)
  if (all(!is.na(variance.delta)))
    checkVectorType(variance.delta)
  if (all(!is.na(innovations.variance)))
    checkVectorType(innovations.variance)
  checkVectorType(delta.range)
  checkVectorType(delta)
  if (!is(dictionary, "wavDictionary"))
    stop("dictionary must be an object of class \"wavDictionary\"")
  checkVectorType(levels, "integer")
  checkScalarType(edof.mode, "integer")
  if (!is.list(boundary) || !all(is.element(names(boundary), c("mode","description"))))
    stop("boundary must be a list with named objects ",
      "\"mode\" and \"description\"")
  checkScalarType(boundary$mode, "logical")
  checkScalarType(boundary$description, "character")
  if (!is(series, "signalSeries"))
    stop("series must be an object of class \"signalSeries\"")
  checkScalarType(sdf.method, "character")
  checkScalarType(type, "character")
  type <- match.arg(type, c("block","instantaneous"))

  z <- list(estimator=estimator,
    delta=delta,
    variance.delta=variance.delta,
    innovations.variance=innovations.variance,
    delta.range=delta.range,
    dictionary=dictionary,
    levels=levels,
    edof.mode=edof.mode,
    boundary=boundary,
    series=series,
    sdf.method=sdf.method,
    type=type)

  oldClass(z) <- "wavFDP"

  return(z)
}

###
# print.wavFDP
##

"print.wavFDP" <- function(x, digits=5, justify="left", sep=":", ...)
{
  # define series string
  name <- ifelse1(!length(x$series@title), "", paste("of", x$series@title))
  type <- lowerCase(x$type)
  est  <- lowerCase(x$estimator)

  if (type == "block"){

    main <- paste("Block-dependent FD parameter estimation", name)

   	# display block dependent results
   	delta     <- x$delta
   	var.delta <- x$variance.delta
   	innov     <- x$innovations.variance

   	# remove certain figures until resolved
   	if (est == "mle")
   	  var.delta <- NA
    else if (est == "wlse")
      innov <- NA

    z <- list(
      "FD parameter estimate (delta)"=round(delta, digits),
      "var{delta} estimate"=round(var.delta, digits),
      "Innovations variance estimate"=round(innov, digits))
  }
  else if (type == "instantaneous"){

    #cat("\nInstantaneous FD parameter estimation", name, "\n\n")

    main <- paste("Instantaneous FD parameter estimation", name)

	  # display block dependent results
	  mean.delta <- mean(x$delta, na.rm=TRUE)
	  mean.var.delta <- mean(x$variance.delta, na.rm=TRUE)
	  mean.innovations <- mean(x$innovations.variance, na.rm=TRUE)

	  # remove certain figures until resolved
	  if (est == "mle")
	    mean.var.delta <- NA
	  else if (est == "lse")
	    mean.innovations <- NA

	  z <- list(
	    "Mean of FD parameter estimates (deltas)"=round(mean.delta, digits),
      "Mean of var{delta} estimates"=round(mean.var.delta, digits),
	    "Mean of innovations variance estimates"=round(mean.innovations, digits))
  }
  else{

  	main <- "Unknown method."
  }

  dict <- x$dictionary

  z <- c(z, list(
    "Estimator"=upperCase(x$estimator),
    "Levels"=x$levels,
    "Boundary mode"=x$boundary$description,
    "EDOF mode"= ifelse1(est == "wlse" && type == "block", x$edof, NULL),
    "Chi-squared DOF"=ifelse1(type == "instantaneous", x$edof * 2 + 1, NULL),
    "Delta range"=x$delta.range,
    "Wavelet"=dict$wavelet,
    "Length of series"=dict$n.sample,
    "Number of levels"=dict$n.levels,
    "Boundary correction rule"=dict$boundary,
    "Filtering technique"=ifelse1(dict$conv,"Convolution","Correlation")))

  prettyPrintList(z, header=main, sep=sep, justify=justify, ...)
}

###
# summary.wavFDP
##

"summary.wavFDP" <- function(object, ...)
{
  obj  <- object
  oldClass(obj) <- "summary.wavFDP"
  return(obj)
}

###
# print.summary.wavFDP
##

"print.summary.wavFDP" <- function(x, ...)
{
  print.wavFDP(x, ...)
}

###
# plot.wavFDP
##

"plot.wavFDP" <- function(x, mean.delta=NULL, xlab="Time",
  ylab=paste(ifelse1(x$boundary$mode, "Biased ", "Unbiased "), upperCase(x$estimator),"(delta)",sep=""),
  title.str=NULL, type="l", show.key=TRUE, conf.color="lightblue", ...)
{
  # save plot parameters and restore upon exit
  old.plt <- splitplot(1,1,1)
  on.exit(par(old.plt))

  # local functions
  # shade the area between two discrete curves
  "patchLimits" <- function(x, y, t=seq(x),  color=13, border=FALSE, ...)
  {
    # perform argument checking
    if (length(x) != length(y) | length(x) != length(t))
      stop("All input vectors must be the same length")
    if (any(is.na(c(x,y,t))))
      stop("NA value not allowed")

    # concatenate the vertices in the proper order
    xdata <- c(t, rev(t))
    ydata <- c(x, rev(y))

    # plot the points without actually plotting the points
    plot(xdata, ydata,type="n", ...)

    # apply polygon shaded area
    polygon(xdata, ydata, col=color, border=border)
  }

  # check for instantaneous
  if (x$type != "instantaneous")
	stop("The plot method for FDP class objects is only supported for instantaneous estimators")

  # obtain the FDP parameters
  delta <- x$delta
  stdev <- sqrt(x$variance.delta);
  bad   <- which (is.na(delta));
  good  <- which (!is.na(delta));

  # define local plot variables
  pch.bad <- 5
  refstr <- ifelse1(is.numeric(mean.delta), paste(", ref=", mean.delta, sep=""), "")

  # for now this routine supports LSE biased/unbiased and MLE unbiased
  use.confidence <- !x$boundary$mode & (lowerCase(x$estimator) == "lse")

  # plot the error bars
  if (use.confidence){

    if (is.numeric(mean.delta)){

      high <- mean.delta + stdev[good]
      low  <- mean.delta - stdev[good]
    }
    else{

      low  <- delta[good] - stdev[good]
      high <- delta[good] + stdev[good]
    }

    ylim.low  <- min(low, delta[good])
    ylim.high <- max(high, delta[good])
    ylim      <- c(ylim.low, ylim.high)

    patchLimits(low, high, t=good, xlab=xlab, ylab=ylab, ylim=ylim, color=conf.color)

    # plot delta
    lines(good, delta[good], type=type,  col=1)

    # plot boundary coefficient points
    points(bad, rep(0, length(bad)), pch=pch.bad, col="red")

    # add a legend
    if (show.key){

    	   legend("bottom",
           c(paste("Delta(t)", refstr, sep=""),"NA Values","95% Confidence          "),
    	     lty=c(1,-1,-1),
    	     pch=c(-1,pch.bad,15),
    	     col=c(1,"red",conf.color),
           cex=1.5)
    	
    }
  }
  else{

    # plot delta
    plot(good, delta[good], type=type,  col=1, xlab=xlab, ylab=ylab)

    # plot boundary coefficient points
    points(bad, rep(0, length(bad)), pch=pch.bad, col="red")

    # add a mean reference line
    if (is.numeric(mean.delta)){

      # plot mean reference line
      lines(range(good), rep(mean.delta, 2), lty= 2,  col=3)

      # add a legend
      if (show.key){

     	   legend("bottom",
           c("Delta(t)","NA Values", paste("Reference delta:", mean.delta,"               ")),
    	     lty=c(1,1,2),
    	     pch=c(-1,pch.bad,-1),
    	     col=c(1,"red",conf.color),
           cex=1.5)
       }
    }
    else{
      if (show.key){

     	   legend("bottom",
           c("Delta(t)","NA Values"),
    	     lty=c(1,-1),
    	     pch=c(-1,pch.bad),
    	     col=c(1,"red"),
           cex=1.5)
      }
    }
  }

  # add a title
  if (is.null(title.str) && (x$type == "instantaneous"))
     title.str <- paste("Instantaneous FD parameter estimation of", wavTitle(x))

  if (!is.null(title.str))
    title(title.str, cex=0.7)

  invisible(NULL)
}

###
# eda.plot.wavFDP
##

"eda.plot.wavFDP" <- function(x, mean.delta=NULL, xlab="Time",
  ylab=paste(ifelse1(x$boundary$mode, "Biased ", "Unbiased "), x$estimator,"(delta)",sep=""),
  title.str=NULL, type="l", ...)
{
  # local functions
  "edatitle" <- function(x) title(paste("\n",x,sep=""), adj=1, cex=0.6)

  # obtain the FDP parameters
  delta <- x$delta
  stdev <- sqrt(x$variance.delta);
  bad   <- which (is.na(delta));
  good  <- which (!is.na(delta));

  if (is.null(mean.delta))
    mean.delta <- mean(delta[good])

  # set up a 2x2 grid of plots
  old.plt <- splitplot(2,2,1)
  on.exit(par(old.plt))

  # create lines for the histogram ...
  delta.std <- unique(stdev)
  yloc      <- c(mean.delta - delta.std, mean.delta + delta.std)

  # plot the delta results with confidence intervals
  plot(x, show.key=FALSE, mean.delta=mean.delta)
  edatitle(paste("Delta(t) with 95% Confidence, ref=", round(mean.delta, 5), sep=""))
  ylim <- par("yaxp")

  # plot the probability density ...
  splitplot(2,2,2)
  hist(delta, prob=TRUE, xlab=ylab, horiz=TRUE,  xpd=TRUE, ylim=ylim, col=2)
  edatitle("Probability Density")
  box(1)

  xrange <- par("usr")[1:2]

  # ... draw mean line
  segments(xrange[1], mean.delta, xrange[2], mean.delta, lwd=3, col=8)

  # ... draw confidence interval limits
  xlocmin <- rep(xrange[1], length(yloc))
  xlocmax <- rep(xrange[2], length(yloc))

  segments(xlocmin, yloc, xlocmax, yloc, lty=2)

  # plot the delta results with confidence intervals
  splitplot(2,2,3)
  plot(x, show.key=FALSE, mean.delta=NULL)
  edatitle("Delta(t) with 95% Confidence")
  ylim <-  par("yaxp")

  # plot qqnorm and qqline
  splitplot(2,2,4)
  ylab <- paste(ifelse1(x$boundary$mode, "Biased ", "Unbiased "), upperCase(x$estimator),"(delta)",sep="")
  qqnorm(delta[good], ylab=ylab, mkh=0.02, ylim=ylim)
  qqline(delta[good], col=8)
  edatitle("QQ Normal")

  invisible(NULL)
}
