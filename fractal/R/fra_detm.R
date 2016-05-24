###################################################
## FRACTAL determinism functionality
##
##::::::::::::::::::::::::::::::::::::::::::::::
##
## Class: determinism
## Constructor function: determinism
## Methods:
##
##   eda.plot.determinism
##   plot.determinism
##   print.determinism
##   summary.determinism
##
###################################################

###
# determinism (constructor)
###

"determinism" <- function(x, dimension=6, tlag=NULL, olag=1,
  scale.min=NULL, scale.max=NULL, resolution=NULL,
  method="ce", n.realization=10, attach.summary=TRUE, seed=0)
{
  # define local functions
  "deltaEpsilon" <- function(x, dimension, tlag, olag,
    scale.min, scale.max, resolution, trim=FALSE){

  	# set the image lag to unity since we rely primarily on the
  	# orbital lag for separation images from the origins
  	ilag <- 1

  	# calculate the delta-epsilon statistics
    z <- itCall( "RS_fractal_determinism_delta_epsilon",
	    x, dimension, tlag, olag, ilag,
	    scale.min, scale.max, resolution, trim)
	    #COPY=rep(FALSE,9),
	    #CLASSES=c("matrix", rep("integer",4), rep("numeric",3), "logical"),
	    #PACKAGE="ifultools")

    scale <- z[[1]]
    e     <- z[[2]]
    E     <- z[[3]]

    bad <- which(E == 0.0)
    if (length(E) > 0)
      E[bad] <- NA

    dimstr   <- paste("E=",seq(dimension),sep="")
    scalestr <- paste("r=",round(scale,3),sep="")
    nms      <- list(scalestr, dimstr)

    dimnames(e) <- nms
    dimnames(E) <- nms

    list(e=e, E=E, scale=scale)
  }

  # obtain series name
  data.name <- deparseText(substitute(x))

  # perform checks on input arguments
  if (!is.numeric(x) && !is(x, "signalSeries"))
    stop ("Input data must be numeric")
  if (dimension < 2)
    stop("Maximal embedding dimension must be greater than two")

  # define defaults for missing inputs
  x   <- as.double(as.numeric(x))
  bad <- which(is.na(x))
  if (length(bad > 0))
    x <- x[-bad]
  dx <- diff(range(x))

  if (is.null(scale.max) || is.null(scale.min)){

    mindx <- min(diff(sort(x)))

    if (is.null(scale.min))
      scale.min <- mindx / 1000;
    if (is.null(scale.max))
      scale.max <- dx * sqrt(dimension)
  }

  if (is.null(resolution))
    resolution <- dx / 1000
  if (is.null(tlag))
    tlag <- timeLag(x, method="acfdecor")

  # check input argument types
  checkScalarType(dimension,"integer")
  checkScalarType(tlag,"integer")
  checkScalarType(olag,"integer")
  checkScalarType(scale.min,"numeric")
  checkScalarType(scale.max,"numeric")
  checkScalarType(resolution,"numeric")
  checkScalarType(method,"character")
  checkScalarType(n.realization,"integer")
  checkScalarType(seed,"integer")

  # check input argument values
  if (olag < 1)
    stop("The orbital lag must be positive")
  if (resolution < 0.0)
    stop("The resolution must be positive")
  if (scale.min < 0.0)
    stop("The minimum scale must be positive")
  if (scale.max < 0.0)
    stop("The maximum scale must be positive")
  if (scale.min >= scale.max + resolution)
    stop("The minimum scale must be greater than the maximum scale plus the one scale increment (resolution)")
  if (n.realization < 1)
    stop("Number of surrogate time series realizations must be positive")

  # calculate the delta-epsilon and E-statistics for the original data
  original <- deltaEpsilon(x,
    dimension=dimension,
    tlag=tlag,
    olag=olag,
    scale.min=scale.min,
    scale.max=scale.max,
    resolution=resolution,
    trim=FALSE)

  # create random seeds
  if (seed > 0){
    set.seed(seed)
    seeds <- sample(.Machine$integer.max/100, size=n.realization)
  }
  else
    seeds <- rep(0,n.realization)

  # calculate the results for the bootstrap ensemble
  # do not trim so that it will be easy to align
  # the results in scale. otherwise different ranges
  # of scale could come back for different surrogates,
  # although in all cases the resolution (scale interval)
  # would remain constant
  bootstrap <- lapply(seq(n.realization),
    function(i, x, dim, tlag, olag, rmin, rmax, dr, method, deltaEpsilon, seeds){

      z <- surrogate(x, method=method, seed=seeds[i])

      deltaEpsilon(as.vector(z), dim=dim,tlag=tlag,olag=olag,
        scale.min=rmin, scale.max=rmax, resolution=dr, trim=FALSE)
    },
    x=x,
    dim=dimension,
    tlag=tlag,
    olag=olag,
    rmin=scale.min,
    rmax=scale.max,
    dr=resolution,
    method=method,
    deltaEpsilon=deltaEpsilon,
    seeds=seeds)

  z <- list(original=original, bootstrap=bootstrap)
  oldClass(z) <- "determinism"

  method <- match.arg(lowerCase(method), c("phase","aaft","dh","ce"))

  attr(z, "data.name")  <- data.name
  attr(z, "n.sample")   <- length(x)
  attr(z, "dimension")  <- dimension
  attr(z, "tlag")       <- tlag
  attr(z, "olag")       <- olag
  attr(z, "scale.min")  <- scale.min
  attr(z, "scale.max")  <- scale.max
  attr(z, "resolution") <- resolution
  attr(z, "n.realization") <- n.realization
  attr(z, "method") <- switch(method,
    phase="Phase Randomization",
    aaft="Theiler's Amplitude Adjusted Fourier Transform",
    dh="Davison-Hinkley",
    ce="Circulant Embedding")
  if (attach.summary)
    attr(z,"summary") <- summary(z)
  z
}

###
# eda.plot.determinism
###

"eda.plot.determinism" <- function(x, summary.=NULL, dimension=NULL, n.scale=NULL, n.box=NULL,
  density=c(0,10), angle=c(45,-45), col=2:3, add=FALSE,
  main=paste("Determinism Test for", attr(x,"data.name"), "series"),
  ylab="Determinism Level (%)", ...)
{
  if (!add){
    old.plt <- splitplot(1,1,1)
    on.exit(par(old.plt))
  }

  # calculate boxplot statistics
  if (!is.null(attr(x,"summary")))
    z <- attr(x,"summary")
  else if (is.null(summary.) || !is(summary.,"summary.determinism"))
    z <- summary(x, n.scale=n.scale, n.box=n.box, dimension=dimension)
  else
    z <- summary.

  nms <- as.vector(rbind(paste("     DIM=",seq(z$dimension),sep=""),rep("",z$dimension)))

  y <- barplot(t(as.matrix(data.frame(z$pass))),
    names=nms,
    density=density,
    angle=angle,
    beside=TRUE,
    space=c(0,2),
    ylim=c(0,115),
    col=col, ylab=ylab)
  if (!is.null(main))
    title(main, adj=0)

  loc <- par("usr")[c(1,4)] * c(1,1)
  legend(loc[1], loc[2], legend=c("Quartile","Extreme"), density=density,
    angle=angle, fill=col)

  invisible(NULL)
}

###
# plot.determinism
###

"plot.determinism" <- function(x, summary.=NULL,
  n.scale=NULL, n.box=NULL, dimension=NULL,
  col=1, range=0, whisklty=1, boxcol=1,
  lwd=2, pch=18, cex=0.8, type="b", ...)
{

  old.plt <- splitplot(1,1,1)
  on.exit(par(old.plt))

  # calculate boxplot statistics
  if (!is.null(attr(x,"summary")))
  {
    z <- attr(x,"summary")
  } else if (is.null(summary.) || !is(summary.,"summary.determinism")) {
    z <- summary(x, n.scale=n.scale, n.box=n.box, dimension=dimension)
  } else {
    z <- summary.
  }

  nc <- floor(sqrt(z$dimension))
  nr <- ceiling(z$dimension/nc)

  for (i in seq(z$dimension)){

  	if (i == 1)
    {
  	  old.plt <- splitplot(nr, nc, i)
  	  on.exit(par(old.plt))
  	} else {
  	  splitplot(nr, nc, i)
    }

    main <- paste(
      ifelse1(z$pass$extreme[i] == z$pass$quartile[i], "", paste(z$pass$extreme[i],"%-", sep="")),
        z$pass$quartile[i], "% PASS", sep="")

    bxp(z$E.bxp[[i]], ylim=z$ylim[[i]],
      col=col, range=range, whisklty=whisklty, boxcol=boxcol,
      ylab=ifelse1((i %% nc) == 1, "E",""),cex=1)

    mtext(paste("DIM", i), side=3, adj=0, line=0.2, cex=1)
    mtext(main, side=3, adj=1, line=0.2, cex=1)
    if (ceiling(i/nc) == nr)
      mtext("SCALE", side=1, line=1, cex=1)

    lines(z$iE, z$E.orig.box[[i]], lwd=lwd, pch=pch, cex=cex, type=type, lty=1)
    axis(side=1, at=z$iE, labels=FALSE)
  }

  # add title to page of plots
  mtext(outer=TRUE, attr(x,"data.name"), cex=1.1, col=2, side=3,adj=0)

  invisible(z)
}

###
# print.determinism
###

"print.determinism" <- function(x, justify="left", sep=":", digits=3, ...)
{
  xatt  <- attributes(x)
  main  <- paste("Determinism test for", xatt$data.name, "series")

  z <- list(
    "Number of surrogate realizations"=xatt$n.realization,
    "Surrogate method"=xatt$method,
    "Embedding dimension(s)"=seq(xatt$dimension),
    "Length of series"=xatt$n.sample,
    "Time lag"=xatt$tlag,
    "Orbital lag"=xatt$olag,
    "Scale range"=paste(round(xatt$scale.min,digits), round(xatt$scale.max,digits), sep=" to "),
    "Scale resolution"=xatt$resolution
    )

  prettyPrintList(z, header=main, sep=sep, justify=justify, ...)

  if (!is.null(attr(x,"summary"))){
    cat("\nDeterminism Level Percentages\n(E-statistic Ensemble vs. Embedding Dimension):\n\n")
  	print(t(as.matrix(data.frame(attr(x,"summary")$pass))))
  }

  invisible(x)
}

###
# summary.determinism
###

"summary.determinism" <- function(object, dimension=NULL, n.scale=NULL, n.box=NULL, ...)
{
  if (!is.null(attr(object,"summary")) && is(attr(object,"summary"),"summary.determinism"))
    return(attr(object,"summary"))

  # define local functions
  "passPercentage" <- function(x, conf){
    z <- x > conf[1,] | x < conf[2,]
    bad <- which(is.na(z))
    if (length(bad))
      z[bad] <- NULL
    round(length(which(z))/(length(x)-length(bad))*100,2)
  }

  # assign defaults
  if (is.null(dimension))
    dimension <- ncol(object$original$E)
  if (is.null(n.scale))
    n.scale <- min(50, length(object$original$scale))
  if (is.null(n.box))
    n.box <- 20

  # check input args
  if (n.box < 1)
    stop("n.box must be positive")
  if (n.scale < 1)
    stop("n.scale must be positive")
  if (n.box > n.scale)
    n.box <- n.scale

  # form pseudo x-axis corresponding to boxplots
  ibox <- floor(seq(1, n.scale, length=n.box))
  fac  <- 100 / (3*n.box*(n.box+1))
  from <- fac*(2*n.box+1)
  by   <- fac*(3*n.box+2)
  iE   <- seq(from=from, by=by, length=n.box)

  # allocate space for output
  ylim <- E.boot.box <- E.boot.all <- E.orig.box <- E.orig.all <- E.bxp <- vector("list", dimension)
  pass <- list(quartile=vector("numeric",dimension),
    extreme=vector("numeric",dimension))

  for (i in seq(dimension)){

    # extract data corresponding to the current dimension.
    # put the data into a list such that each object of the
    # list contains a vector corresponding to the values of
    # all the bootstrapped E(dimension,scale) statistics ...

    # form n.scale by n.realization matrix containing the
    # surrogate E values for the current embedding dimension
    E.boot <- sapply(object$bootstrap, function(x, i) x$E[,i], i=i)

    # obtain original E-statistic at current dimension
    E.orig <- object$original$E[,i]

    # find common minimum scale at which there are no NA values
    # for the bootstrapped and original E-statistics
    imin.boot <- min(which(apply(E.boot,MARGIN=1,function(x) !any(is.na(x)))))
    imin.orig <- min(which(!is.na(E.orig)))
    imin      <- max(imin.boot, imin.orig)

    # starting from the commom minimum usable scale,
    # form scale index vector
    ibox  <- seq(imin, imin + n.scale, length=n.box)
    iall  <- seq(imin, length=n.scale)
    scale <- object$original$scale[ibox]

    # sample the appropriate rows of the
    # E.boot(n.scale, n.realization) matrix
    # to collect bootstrapped E statisitics for
    # each boxplot. do the same for the original data.
    E.boot.box[[i]] <- lapply(ibox, function(i, m) m[i,], m=E.boot)
    E.boot.all[[i]] <- lapply(iall, function(i, m) m[i,], m=E.boot)
    E.orig.box[[i]] <- E.orig[ibox]
    E.orig.all[[i]] <- E.orig[iall]

    # define y-axis limits to avoid plot warnings
    yrange <- range(unlist(E.boot.box[[i]]), E.orig.box[[i]], na.rm=TRUE)
    dy     <- 0.05 * abs(diff(yrange)) * c(-1, 1)
    ylim[[i]] <- yrange + dy

    # collect boxplot statistics for ALL
    # data within usable range of scales
    boxstats <- boxplot(E.boot.all[[i]],plot=FALSE)
    conf     <- boxstats$conf
    extr     <- boxstats$stats[c(1,5),]

    # calculate pass-fail summaries for the E-statistics
    # in the current dimension, one for comparison to
    # the first quartiles of the bootstrapped statistics,
    # and one for the extremes (excluding outliers)
    pass$quartile[i] <- passPercentage(E.orig.all[[i]], conf)
    pass$extreme[i]  <- passPercentage(E.orig.all[[i]], extr)

    # create scale labels for the x-axis and main title
    scale.names <- rep("", n.box)
    scale.names[1] <- as.character(round(scale[1],3))
    scale.names[n.box] <- as.character(round(scale[length(scale)],3))

    # create the boxplot data
    E.bxp[[i]] <- boxplot(E.boot.box[[i]], names=scale.names, plot=FALSE)
  }

  z <- list(E.boot.box=E.boot.box, E.boot.all=E.boot.all, E.orig.box=E.orig.box,
    E.orig.all=E.orig.all, E.bxp=E.bxp, iE=iE, ylim=ylim, pass=pass,
    n.box=n.box, n.scale=n.scale, dimension=dimension)

  oldClass(z) <- "summary.determinism"
  z
}

