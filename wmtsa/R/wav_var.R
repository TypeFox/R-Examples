################################################
## WMTSA package wavelet variance functionality
##
##  Functions:
##
##    D.statistic
##    D.table
##    wavVarConfidence
##    wavEDOF
##
##  Constructor Functions and methods:
##
##    wavVar
##
##      plot.wavVar
##      print.wavVar
##      summary.wavVar
##      print.summary.wavVar
##
##    wavVarTest:
##
##      plot.wavVarTest
##      print.wavVarTest
##
##################################################

###
# D.statistic
##

"D.statistic" <- function(x)
{
  # The D.statistic() function calculates the
  # the maximum departure from an expected
  # linear increase in cumulative energy

  N      <- length(x)
  Nm1    <- (N - 1)
  cumulative.energy <- cumsum(x[1:(length(x)-1)]^2)
  P      <- cumulative.energy / sum(x^2)
  index  <- seq(length=Nm1, from=0)
  Dplus  <- max((index + 1) / Nm1 - P)
  Dminus <- max(P - (index / Nm1))
  D      <- max(c(Dplus, Dminus))
  # my (dbp) very first comment!!! This was a great comment!
  D
}

###
# D.table
##

"D.table" <- function(n.sample=c(127, 130), significance=c(0.1,0.05,0.01),
  lookup=TRUE, n.realization=10000, n.repetition=3, tolerance=1e-6)
{
  # define local functions
  "D.table.preprocess" <- function(significance, n.sample)
  {
    # The D.table.preprocess() function sorts the significance
    # and n.sample vectors and checks for appropriate ranges.
    # The new n.sample and significance vectors are returned
    # in a list object.

    # test to make sure the significance is in the proper range
    significance.range <- range(significance)
    checkRange(significance.range)

    n.sample <- as.integer(n.sample)

    if (!(all(n.sample > 0)))
      stop("n.sample value(s) must be positive integers")

    # sort the significances from high to low
    significance <- unique(rev(sort(significance)))

    # sort the n.sample vector from low to high
    n.sample <- unique(sort(n.sample))

    list(significance=significance, n.sample=n.sample)
  }
 
  assign("D.table.preprocess", D.table.preprocess)
 
  "D.table.names" <- function(n.sample, significance)
  {
    # The D.table.names() function creates the names for the
    # rows and columns of a D.table matrix
    ilow  <- which(n.sample < 128)
    method <- rep("Inclan-Tiao",length(n.sample))
    method[ilow] <- rep("Monte Carlo",length(ilow))
    list(method , c("N", paste(significance*100,"%",sep="")))
  }

  assign("D.table.names", D.table.names)

  "significance.percentage.to.number" <- function(x)
  {
    # The significance.percentage.to.number() function
    # converts the significance(s) expressed as a
    # character vector with a percentage to a number(s)
    # For example, "10%" goes to 10.

    # get rid of "%" characters and turn column
    # names into significances
    as.numeric(gsub("%","",x))/100
  }

  "D.table.make" <- function(significance=c(0.1, 0.05, 0.01), n.sample=c(10,200),
    n.repetition=3, n.realization=10000, tolerance=1e-6)
  {
    # The D.table.make() function returns a matrix of D-statistics
    # corresponding to a single round of Monte Carlo simulations
    # if n.sample < 128 combined with the D-statistics
    # calculated via Incan-Tiao for n.sample >= 128.
    # The rows are the realization lengths while the columns are
    # the significances.

    # define local functions

    "D.MonteCarlo" <- function(n.sample=10, n.realization=10000, significance=0.1)
    {
      # The D.MonteCarlo() function returns the D-statistics corresponding to
      # the significance(s) for a single realization length (n.sample < 128)

      # precalculate the indices of the sorted D-statistic vector
      # corresponding to the significances
      index <- floor((1 - significance) * n.realization)

      D.realization <- sort(unlist(lapply(seq(length=n.realization),
        function(i, n.sample){

          # generate a N(0,1) realization
          G <- rnorm(n=n.sample, mean=0, sd=1)

          # calculate and return the D-statistic
          D.statistic(G)
        },
        n.sample=n.sample)))

      # return D-statistics corresponding to the significances
      (D.realization[index] + D.realization[index + 1]) / 2
    }

    "D.InclanTiao" <- function(n.sample=10, tolerance=1e-6, significance=0.1)
    {
      # define local functions
      "InclanTiao.root" <- function(d, n.sample=100, tolerance=1e-6, significance=0.1)
      {
        # The Inclan-Tiao root function is used in
        # combination with the uniroot function to approximate
        # the D-statistic for sample sizes larger than (say) 128.
        # The true function involves an infinite summation of exponentials
        # which tend to decay rather quickly. Thus, a tolerance
        # denoting the amplitude threshold limits the number of
        # of summations and increases the speed of calculation.

        if (n.sample > 100){
          max.iterations <- ceiling(sqrt(-log(tolerance/2) / n.sample / d^2))
          l <- seq(length=max.iterations)
          2 * sum((-1)^(l+1) * exp(-l*l*n.sample*d*d)) - significance
        }
        else
          stop("The Inclan-Tiao method works best for N.SAMPLE > 100")
      }

      # returns the D-statistics corresponding to
      # the significance(s) for a single realization length (n.sample >= 128)

      unlist(lapply(seq(along=significance),
	    function(i, n.sample, significance, tolerance,InclanTiao.root){
          eps <- .Machine$double.eps
          uniroot(InclanTiao.root,
            c(1e-4,1-eps), n.sample=n.sample,
            significance=significance[i],
            tolerance=tolerance)$root
          },
          n.sample=n.sample,
		  significance=significance,
		  tolerance=tolerance,
		  InclanTiao.root=InclanTiao.root))
    }

    # PRIMARY CODE SET
    # check ranges on input and sort
    input        <- D.table.preprocess(significance, n.sample)
    significance <- input$significance
    n.sample     <- input$n.sample

    D.nsample <- unlist(lapply(seq(along=n.sample),
      function(i, n.sample, n.realization, significance, tolerance,
        n.repetition, D.MonteCarlo, D.InclanTiao){

        # calculate the critical levels corresponding to
        # each realization length (n.sample)
        if (n.sample[i] < 128){

          D.round <- lapply(seq(length=n.repetition),
            function(j, n.sample, n.realization, significance, D.MonteCarlo){

              D.MonteCarlo(n.sample=n.sample,
                n.realization=n.realization,
                significance=significance)
            },
            n.sample=n.sample[i],
            n.realization=n.realization,
            significance=significance,
            D.MonteCarlo=D.MonteCarlo)

          D.array <- array(unlist(D.round), c(length(significance), n.repetition))
          D <- rowMeans(D.array)
        }
        else{

          D <- D.InclanTiao(n.sample=n.sample[i],
            tolerance=tolerance,
            significance=significance)
        }

        return(D)
      },
      n.sample=n.sample,
      n.realization=n.realization,
      significance=significance,
      tolerance=tolerance,
      n.repetition=n.repetition,
      D.MonteCarlo=D.MonteCarlo,
      D.InclanTiao=D.InclanTiao))

    D.tab <- cbind(n.sample, matrix(D.nsample, ncol=length(significance), byrow=TRUE))

    dimnames(D.tab) <- D.table.names(n.sample, significance)

    return(D.tab)
  }

  "D.table.fill.missing" <- function(D.row, D.table.make, col.names=NULL)
  {
    # The D.table.fill.missing() function fills in the
    # the D-statistics missing from a D.table row

    # define local functions
    "significance.percentage.to.number" <- function(x)
    {
      # The significance.percentage.to.number() function
      # converts the significance(s) expressed as a
      # character vector with a percentage to a number(s)
      # For example, "10%" goes to 10.

      # get rid of "%" characters
      x <- gsub("%","",x)

      # ... turn them into significances
      as.numeric(x)/100
    }

    missing.indices <- is.na(D.row)

    if (any(missing.indices)){

      significance.na <- significance.percentage.to.number(col.names[missing.indices])
      n.sample.na <- D.row[1]
      D.na <- D.table.make(significance=significance.na,
        n.sample=n.sample.na,
        n.realization=100)

      # invoke column offset for significances returned by the
      # D.table.make() function (the first column denotes the n.sample
      # value, so we skip this column when combining results)
      D.na.indices <- seq(from=2,length=length(which(missing.indices)))
      D.row[missing.indices] <- D.na[D.na.indices]
    }

    D.row
  }

  if (lookup){

    # check ranges on input and sort
    input        <- D.table.preprocess(significance, n.sample)
    significance <- input$significance
    n.sample     <- input$n.sample

    # verify existence of critical D statistic table.
    # if not available, issue error
    if (!exists("D.table.critical")){

      warning("The D.table.critical lookup table does not exist. The table is calculated explicitly")

      D <- D.table.make(significance=significance,
        n.sample=n.sample,
        n.repetition=n.repetition,
        n.realization=n.realization,
        tolerance=tolerance)
    }
    else{

      col.names <- dimnames(D.table.critical)[[2]]

      # find column names that have a "%" in them
      ipercentage <- grep("%", col.names)
      col.names   <- col.names[ipercentage]

      # ... turn them into significances
      table.significance <- significance.percentage.to.number(col.names)

      # find the column indices in D.critical which match the
      # significances. Non-matching elements of the significance
      # vector are returned as NA's
      column.offset <- min(ipercentage) - 1
      cols <- match(significance, table.significance) + column.offset

      # find the row indices in in D.critical which match the
      # n.sample values. Non-matching elements of the n.sample
      # vector are returned as NA's
      rows <- match(n.sample,  D.table.critical[,1])

      # build critical D-statistic table
      D <- cbind(N=n.sample, D.table.critical[rows,cols, drop=FALSE])
      dimnames(D) <- D.table.names(n.sample, significance)

      # if there are missing values in the table, fill them in
      if (any(is.na(D))){

	D <- t(apply(D, MARGIN=1,
	  function(x, D.table.fill.missing, D.table.make, col.names){
	    D.table.fill.missing(x, D.table.make, col.names)
	  },
	  D.table.fill.missing=D.table.fill.missing,
	  D.table.make=D.table.make,
	  col.names=dimnames(D)[[2]]))
      }

    } # end check for existence of D.table.critical

  } # end if lookup
  else{

    D <- D.table.make(significance=significance,
		 n.sample=n.sample,
		 n.repetition=n.repetition,
		 n.realization=n.realization,
		 tolerance=tolerance)
  }

  dimnames(D) <- D.table.names(n.sample, significance)

  # remove local functions
  remove("D.table.preprocess")
  remove("D.table.names")

  return(D)
}

###
# wavVarConfidence
###

"wavVarConfidence" <- function(wvar, edof, probability=0.95)
{

  # turn EDOF list into a vector, one EDOF element per scale

  edof <- unlist(edof)

  if (any(is.na(edof))){
    na.conf    <- rep(NA, length(edof))
    confidence <- list(low=na.conf, high=na.conf)
  }
  else{

    confidence <- itCall("RS_wavelets_variance_confidence",
      wvar,
      edof,
      probability)
    #
      #COPY=rep(FALSE,3),
      #CLASSES=c("numeric","numeric","numeric"),
      #PACKAGE="ifultools")

    confidence <- lapply(confidence,function(x,crystals){
	y <- as.vector(x)
	names(y) <- crystals
	return(y)
	},
	crystals=names(edof))

    names(confidence) <- c("low","high")
  }

  return(confidence)
}

###
# wavEDOF
###

"wavEDOF" <- function(x, wavelet="s8", levels=NULL, sdf=NULL, sdfargs=NULL, sampling.interval=1, n.fft=1024)
{
  if (is(x,"wavTransform")){

    # obtain interior wavelet coefficients
    Wb <- wavBoundary(x)
    J  <- x$dictionary$n.levels

    # form levels vector if empty
    if (is.null(levels))
      levels <- seq(J)

    if (max(levels) > J)
      stop("Levels exceeds maximum in transform")

    crystals <- paste("d", levels, sep="")
    interior <- Wb$interior[crystals]
    n.coeff  <- Wb$interior$length[crystals]

    # form the unbiased wavelet variance estimates
    variance <- unlist(lapply(interior,function(x) sum(x^2))) / n.coeff

    # adjust normalization for DWT style
    if (x$xform == "dwt")
      variance <- variance / 2^(levels)

    # map filters
    filter <- mutilsFilterType(x$dictionary$wavelet)

    # evaluate input SDF over frequencies [0, Nyquist]
    Sx <- mutilsSDF(sdf, sdfargs=sdfargs, n.freq=n.fft, sampling.interval=sampling.interval)

    # adjust for case where there are no interior coefficients
    # at a given level
    bad <- which(is.na(variance))
    if (length(bad)){
    	interior <- interior[-bad]
    	n.coeff  <- n.coeff[-bad]
    	levels   <- levels[levels != bad]
    	crystals <- crystals[-bad]
    	variance <- variance[-bad]
    }

    # calculate the EDOF
    edof <- itCall("RS_wavelets_variance_edof",
      as.numeric(unlist(interior)),
      as.integer(n.coeff),
      as.numeric(variance),
      as.integer(levels),
      as.numeric(Sx),
      as.integer(filter$type),
      as.integer(filter$length))
    #
      #COPY=rep(FALSE,7),
      #CLASSES=c("numeric","integer","numeric","integer","numeric","integer","integer"),
      #PACKAGE="ifultools")

    edof <- lapply(edof,
      function(x,crystals){
        y <- as.vector(x)
        names(y) <- crystals
        invisible(y)}, crystals=crystals)

    names(edof) <- paste("EDOF", 1:3, sep="")
    names(n.coeff) <- crystals

    # pack results into a list
   	edof <- c(edof, list(variance.unbiased=variance, n.coeff=n.coeff))

    # if any members of edof are < 0, that means that they
    # were not calculated. in this case, place an NA for
    # each element of that edof mode
    if (any (edof[[2]] < 0))
      edof[[2]] <- rep(NA,length(edof[[2]]))

    return(edof)
  }
  else if (is.element(class(x), c("numeric","ts","rts","signalSeries"))){

    y <- create.signalSeries(x)
    dt <- y@positions@by

    if (is.null(levels))
      n.level <- ilogb(length(y),base=2)
    else
      n.level <- as.integer(max(levels))

    xform <- wavMODWT(y, wavelet=wavelet, n.levels=n.level)
    return(wavEDOF(xform, sdf=sdf, levels=levels, sampling.interval=dt, sdfargs=sdfargs))
  }
  else{
    stop("Class of input is unsupported")
  }
}

##############################################################################
##
##  Constructor wavVar:
##
##    plot.wavVar
##    print.wavVar
##    summary.wavVar
##    print.summary.wavVar
##
##############################################################################

###
# Constructor: wavVar
###

"wavVar" <- function(x, xform="modwt", wavelet="s8", n.levels=NULL,
  position=list(from=1,by=1,units=character()), units=character(),
  documentation=character(), sdf=NULL, sdfargs=NULL, sampling.interval=1, n.fft=1024)
{
  # get the series name from the parent
  # (or from this function if this is called as a top level function)
  series.name <- deparse(substitute(x))
  series <- create.signalSeries(x, position=position, units=units,
  title.data=series.name, documentation=documentation)

  # map filters
  filter <- mutilsFilterType(wavelet)

  if (is.null(n.levels))
    n.levels <- wavMaxLevel(n.taps=filter$length, n.sample=length(x), xform=xform)

  # evaluate input SDF over frequencies [0, Nyquist]
  Sx <- mutilsSDF(sdf=sdf, sdfargs=sdfargs, n.freq=n.fft, sampling.interval=deltat(x))

  type <- mutilsTransformType(xform)

  obj <- itCall("RS_wavelets_variance",
    as.numeric(series@data),
    as.integer(type),
    as.integer(filter$type),
    as.integer(filter$length),
    as.integer(n.levels),
    as.numeric(Sx))
    #
    #COPY=rep(FALSE,6),
    #CLASSES=c("numeric","integer","integer","integer","integer","numeric"),
    #PACKAGE="ifultools")

  is.conf <- as.logical(lowerCase(xform) == "modwt")

  # develop crystal names
  levels   <- seq(along=obj[[2]][[1]])
  crystals <- paste("d",levels,sep="")
  sampling.interval <- series@positions@by
  scales   <- 2^(levels - 1)*sampling.interval

  # name the data
  names(scales) <- crystals

  vtime <- lapply(obj[[1]], as.vector)
  names(vtime) <- c("biased","unbiased")

  vblock <- lapply(obj[[2]], function(x,crystals){
	y <- as.vector(x);names(y) <- crystals;return(y) }, crystals=crystals)
  names(vblock) <- c("biased","unbiased")

  filters <- wavDaubechies(wavelet=wavelet, normalized=(xform == "modwt"))

  if (is.conf){

    conf <- lapply(obj[[3]], function(x,crystals){
      low  <- x[1,]
      high <- x[2,]

      names(low)  <- crystals
      names(high) <- crystals

      return(list(low=low,high=high))
    },
    crystals=crystals)

    names(conf) <- paste("n",1:3,sep="")

    separation <- itCall("RS_wavelets_transform_coefficient_boundaries",
      as.integer(n.levels), as.integer(filters$length), as.integer(length(series)), as.integer(type))
    #
      #COPY=rep(FALSE,4),
      #CLASSES=rep("integer",4),
      #PACKAGE="ifultools")

    unbiased.length <- separation[[3]][levels]
    names(unbiased.length) <- crystals

    edof <- lapply(obj[[4]],
      function(x,crystals){
        y <- as.vector(x)
        names(y) <- crystals
        return(y)
        }, crystals=crystals)
    names(edof) <- paste("EDOF",1:3,sep="")
  }
  else{
     conf <- list(EDOF1=NA, EDOF2=NA, EDOF3=NA)
     edof <- list("n1"=NA, "n2"=NA, "n3"=NA)
     unbiased.length <- NA
  }

  # create dictionary
  dict <- list(wavelet =wavelet,
    dual     = FALSE,
    decimate = FALSE,
    n.sample = length(x),
    attr.x   = NULL,
    n.levels = n.levels,
    boundary = "periodic",
    conv     = TRUE,
    analysis.filter=list(low=filters$scaling, high=filters$wavelet),
    synthesis.filter=list(low=filters$scaling, high=filters$wavelet),
    fast=TRUE,
    is.complex=FALSE)

  oldClass(dict) <- "wavDictionary"

  series@title <- series.name

  wvar <- list(time=vtime,
    block=vblock,
    scales=scales,
    xform=xform,
    dictionary=dict,
    series=series,
    sampling.interval=sampling.interval)

  z <- c(wvar, list(confidence=c(edof,conf,list(length=unbiased.length))))

  if (is.null(sdf))
    z$confidence$EDOF2 <- NA

  oldClass(z) <- "wavVar"

  z
}

###
# plot.wavVar
###

"plot.wavVar" <- function(x, type="unbiased", time=FALSE, xlab=NULL,
	title=NULL, ylab=NULL, edof=3, pch=18, logfunc=log10, ...)
{
  # check input data
  checkVectorType(edof,"integer")
  checkRange(edof,c(1,3))

  # now that the internal functions are developed for plotting,
  # begin plot method code ...

  # map the plot type
  biased <- type == "biased"

  # obtain sequence name and, if too long, truncate
  series            <- x$series
  series.name       <- series@title
  sampling.interval <- x$sampling.interval

  if (length(series.name) < 1)
    series.name <- "x"
  if (length(sampling.interval) < 1)
    sampling.interval <- 1

  # obtain transform type, and wavelet type
  xform          <- upperCase(x$xform)
  wavelet        <- x$dictionary$wavelet
  wavelet.str    <- paste("(", wavelet, ")")
  biased.title   <- paste("Biased", xform, "Wavelet Variance", wavelet.str)
  unbiased.title <- paste("Unbiased", xform, "Wavelet Variance", wavelet.str)

  # create matrix of variance data
  biased.data   <- x$block$biased
  unbiased.data <- x$block$unbiased

  # obtain the crystals for each type of data
  unbiased.crystals <- names(unbiased.data)
  biased.crystals   <- names(biased.data)

  if (is.null(xlab)){
    if (time == FALSE) {

      xlab   <- "Scale"
      xunits <- x$series@units.position
      if (length(xunits) > 0)
        xlab <- paste(xlab, " (", xunits, ")",sep="")
    }
    else
      xlab <- "Time"
  }

  if (is.null(title))
    title <- paste("Wavelet Variance of", series.name)

  if (time){

    # time-dependent wavelet variance plots
    old.plt <- splitplot(1,1,1)
    on.exit(par(old.plt))

    same.scale <- TRUE   # make sure magnitudes are similar
    zero.line  <- TRUE

    N     <- x$dictionary$n.sample
    dt    <- x$series@positions@by
    times <- (1:N) * dt

    nms   <- rev(names(x$time$unbiased))

    if (biased){
      wavStackPlot(x$time$biased[nms], times=times, same.scale=same.scale, zeroline=zero.line)
      title(paste(biased.title,"\n\n",series.name),xlab=xlab)
    }
    else{
      wavStackPlot(x$time$unbiased[nms], times=times, same.scale=same.scale, zeroline=zero.line)
      title(paste(unbiased.title,"\n\n",series.name),xlab=xlab)
    }

  }
  else{

  	# block-dependent wavelet variance plot
    xdata <- logfunc(x$scales)
    ydata <- logfunc(ifelse1(biased, biased.data, unbiased.data))

    if (any(is.na(x$confidence$EDOF2)))
      edof <- edof[edof != 2]

    if (!length(edof))
      stop("No data exists for current EDOF specification")

    conf <- ifelse1(xform == "DWT", list(), x$confidence[paste("n", edof, sep="")])

    # use the same plot limits for all EDOF plots for comparison purposes
    xlim <- range(xdata, na.rm=TRUE)
    if (length(conf))
      ylim <- range(logfunc(unlist(conf)), na.rm=TRUE)

    logtext <- deparse(substitute(logfunc))

    xlab <- paste(logtext, "(scale)", sep="")
    ylab <- paste(logtext, ifelse1(biased, biased.title, unbiased.title))

    for (i in seq(length(edof))){

    	old.plt <- splitplot(1,length(edof),i, gap=0.1)
      if (i == 1)
        on.exit(par(old.plt))

    	if (length(conf)){

    	  low  <- logfunc(conf[[i]]$low)
    	  high <- logfunc(conf[[i]]$high)

    	  plot(xlim, ylim, type="n", xlab=xlab, ylab=ifelse1(i==1,ylab,""))
    	  
       	  arrows(xdata, low, xdata, high, code=3, col="black", angle=90, lwd=2, length=0.1)
       	  points(xdata, ydata, pch=19, col="blue", cex=1.5)
 
        mtext(paste("EDOF", edof[i]), side=3, adj=1, line=0.5)
    	}
    	else
    	  plot(xdata, ydata, pch=19, col="blue", cex=1.5)

    	
    	grid()
     }

  }

  invisible(NULL)
}

###
# print.wavVar
###

"print.wavVar" <- function(x, justify="left", sep=":", ...)
{
  series            <- x$series
  title             <- series@title
  sampling.interval <- x$sampling.interval
  pos.units         <- series@units.position

  if (length(title) < 1)
    title <- "x"

  crystals  <- wavSortCrystals(names(x$block$biased))
  dict      <- x$dict
  main      <- paste(upperCase(x$xform)," wavelet variance of ", title, sep="")

  z <- c(as.list(x$dict), list("Sampling interval"=sampling.interval,
    "Units"=ifelse1(length(pos.units) > 1,  paste("(", pos.units, ")"), NULL)))
  prettyPrintList(as.list(x$dict), header=main, sep=sep, justify=justify, ...)

  scales <- x$scales[crystals]
  print(matrix(scales,nrow=1,dimnames=list("Scales",crystals)), digits=2)

  invisible(x)
}

###
# summary.wavVar
###

"summary.wavVar" <- function(object, ...)
{
  # create matrix of variance data
  biased.data   <- object$block$biased
  unbiased.data <- object$block$unbiased

  # obtain the crystals for each type of data
  unbiased.crystals <- wavSortCrystals(names(unbiased.data))
  biased.crystals   <- wavSortCrystals(names(biased.data  ))

  # bind the variance data into a row major matrix.
  # Since the unbiased estimates use only interior
  # wavelet coefficients, and since some crystals
  # may not contain any interior coefficients,
  # those crystals are ignored. Consequently, the
  # length of the biased and unbiased wavelet vectors
  # may differ.
  data <- matrix(NA,nrow=2,ncol=length(biased.data))
  dimnames(data) <- list(c("Biased","Unbiased"),biased.crystals)
  data[1,] <- biased.data[biased.crystals]
  data[2,][unbiased.crystals] <- unbiased.data[unbiased.crystals]

  # obtain sequence name and, if too long, truncate
  series <- object$series
  series.name <- series@title
  sampling.interval <- object$sampling.interval
  if (length(series.name) < 1)
    series.name <- "x"
  if (length(sampling.interval) < 1)
    sampling.interval <- 1

  cut <- 15
  if (nchar(series.name) > cut)
	  series.name <- paste(substring(series.name,1,cut), " ... ",sep="")

  J <- object$dictionary$n.levels

  if (object$xform == "modwt"){

    confidence <- object$confidence

    # produce matrix to hold data
    scale.data  <- unlist(object$scales[unbiased.crystals])
    EDOF1.data  <- unlist(round(confidence$EDOF1[unbiased.crystals]))
    EDOF2.data  <- unlist(round(confidence$EDOF2[unbiased.crystals]))
    EDOF3.data  <- unlist(round(confidence$EDOF3[unbiased.crystals]))
    length.data <- unlist(confidence$length[unbiased.crystals])
    level.data  <- ilogb(scale.data/sampling.interval, base=2) + 1

    conf.mat <- rbind(level.data, scale.data,
		      EDOF1.data, EDOF2.data,EDOF3.data,
		      length.data)

    # supply names for each dimension of matrix
    dimnames(conf.mat) <- list(c("Level","Scale","[EDOF1]","[EDOF2]","[EDOF3]","# coeffs"),
      biased.crystals[1:ncol(conf.mat)])
  }
  else
    conf.mat <- NA

  obj <- list(conf.mat=conf.mat,
		xform=object$xform,
		vmat=data,
		series.name=series.name)

  oldClass(obj) <- "summary.wavVar"
  return(obj)
}

###
# print.summary.wavVar
###

"print.summary.wavVar" <- function(x, digits=max(2, .Options$digits - 4), ...)
{
  main <- paste(upperCase(x$xform)," wavelet variance of ", x$series.name,sep="")
  cat(main,"\n", rep("-",nchar(main)),"\n",sep="")

  print(x$vmat, digits=digits)

  if (x$xform == "modwt"){
    cat("\nConfidence Interval Data for Unbiased\nTime Independent Wavelet Variance:\n\n")
    print(x$conf.mat)
  }
  invisible(x)
}

##############################################################################
##
##  Constructor wavVarTest:
##
##    plot.wavVarTest
##    print.wavVarTest
##
##############################################################################

###
# wavVarTest (constructor)
##

"wavVarTest" <- function(x, wavelet="s8", n.levels=NULL,
  significance=c(0.1,0.05,0.01), lookup=TRUE, n.realization=10000,
  n.repetition=3, tolerance=1e-6)
{
  # define local functions
  "is.dyadic" <- function(x) as.logical(2^ilogb(x,base=2) == x && x > 0)

  "D.wavebound" <- function(x, good.crystals)
  {
  	# Calculate the D statistic, which estimates the largest departure from an assumed
    # linear increase in cumulative (interior wavelet coefficient)
    # energy. If the wavelet energy is homogeneously distributed across time,
    # its cumulative energy should increase approximately linearly.

    if (!is(x,"wavBoundary"))
      stop("Input must be of class \"wavBoundary\"")

    # obtain energy interior coefficients
    # in each crystal
    energy   <- x$interior$energy[good.crystals]

    # obtain nonzero energy components and corresponding
    # crystal names
    energy   <- energy[energy > 0]
    crystals <- names(energy)

    # return D-statistic
    unlist(lapply(x$interior[crystals], D.statistic))
  }

  if (is(x,"wavBoundary")){

    # check for convolution style filtering
    if (!x$dictionary$conv)
      stop("Convolution style filtering must be used in the DWT ",
        "for for homogeneity of variance tests")

    series.name <- wavTitle(x)
  }
  else{

    # get the series name
    series.name <- deparseText(substitute(x))

    # check for power of 2
    dyad <- is.dyadic(length(x))

    # perform a dwt on the data
    n.sample <- length(x)

    if (is.null(n.levels))
	  n.levels <- wavMaxLevel(
        n.taps= mutilsFilterType(wavelet)$length,
        n.sample=n.sample,
        xform="dwt")

    # perform the dwt using convolution style filtering
    W <- wavDWT(x, wavelet=wavelet, n.levels=n.levels)

    # perform a wavelet boundary separation
    WW <- wavBoundary(W)
  }

  # obtain dwt specifications
  dict       <- WW$dictionary
  n.levels   <- dict$n.levels
  series     <- WW$series
  wave.class <- class(W)
  crystals   <- WW$crystals

  # estimate departure from an assumed linear increase in cumulative
  # (interior wavelet coefficient) energy
  Nj <- WW$interior$length
  Nj <- Nj[Nj > 0]
  good.crystals <- names(Nj)

  # get rid of scaling coefficient crystals
  bad <- charmatch("s", good.crystals)

  if (!is.na(bad)){
    Nj <- Nj[-bad]
    good.crystals <- good.crystals[-bad]
  }

  Dj <- D.wavebound(WW, good.crystals)

  # create a D-statistics table
  D.critical <- D.table(n.sample=Nj,
    significance=significance,
    lookup=lookup,
    tolerance=tolerance,
    n.repetition=n.repetition,
    n.realization=n.realization)

  # if duplicate Nj exist, replicate that row in D.critical.
  # this will occur with sJ and dJ wavelet crystals for example
  D.critical <- D.critical[pmatch(Nj,D.critical[,1],duplicates.ok=TRUE),]

  # perform the homogeneity test
  is.homogeneous <- apply(D.critical[,-1], 2,
    function(x,Dj) Dj < x ,Dj=Dj[good.crystals])

  dimnames(is.homogeneous)[[2]] <-
    dimnames(D.critical)[[2]][-1]

  crystal <- names(Dj)

  # explicitly set the series name
  series@title <- series.name

  z <- list(pass=is.homogeneous,
    D=Dj,
    D.critical=D.critical,
    significance=significance,
    n.sample=Nj,
    dictionary=dict,
    series=series,
    lookup=lookup,
    tolerance=tolerance,
    n.repetition=n.repetition,
    n.realization=n.realization)

  oldClass(z) <- "wavVarTest"

  return(z)
}

###
# plot.wavVarTest
##

"plot.wavVarTest" <- function(x, ...)
{
  series.name <- wavTitle(x)
  sigcols <- 2:ncol(x$D.critical)
  key.names <- dimnames(x$D.critical)[[2]][sigcols]

  xrange <- range(x$D.critical[,1])
  yrange <- range(as.vector(x$D.critical[,2:ncol(x$D.critical)]))

  for (i in 1:sigcols[length(sigcols)]){
    if (i == 1)
      plot(x$n.sample, x$D.critical[,i], ylab="D-statistic",
        xlab="Number of Samples", xlim=xrange,ylim=yrange,log="x")
    else
      points(x$n.sample, x$D.critical[,i], pch=i)
  }

  # obtain axis limits [x1 x2 y1 y2]
  ax <- par("usr")

  # add a title
  title(series.name)
}

###
# print.wavVarTest
##

"print.wavVarTest" <- function(x, ...)
{
  series.name <- wavTitle(x)
  dict        <- x$dictionary
  cat("Homogeneity Test for Discrete Wavelet Transform of", series.name, "\n")
  cat("\nPass:\n")
  print(x$pass)
  cat("\nD-statistic critical values comparison:\n")
  methods <- dimnames(x$D.critical)[[1]]

  D.cat <- cbind(N=x$D.critical[,1], "D"=x$D, x$D.critical[,-1])
  dimnames(D.cat)[[1]] <- paste(dimnames(x$pass)[[1]], methods, sep=": ")
  print(D.cat)
  cat("\nLookup:", x$lookup,"\n")

  if (any(dimnames(x$D.critical)[[1]] == "Inclan-Tiao"))
    cat("Inclan-Tiao tolerance: ", x$tolerance,"\n")

  if (x$lookup){
    cat("Repetitions: 3 (lookup),", x$n.repetition, "(non-lookup)\n")
    cat("Realizations/repetition: 10000 (lookup),", x$n.realization, "(non-lookup)\n")
  }
  else{
    cat("Repetitions:", x$n.repetition, "\n")
    cat("Realizations/repetition:", x$n.realization, "\n")
  }

  cat("\n")
  print(dict)
}

