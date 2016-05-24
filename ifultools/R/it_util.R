################################################
# IFULTOOLS package utility functions
#
# ACVStoPACS
# aggregateData
# checkRange
# checkScalarType
# checkVectorType
# decibel
# em
# ilogb
# isVectorAtomic
# linearFit
# linearSegmentation
# logScale
# mergeList
# nextDyadic
# ordinal
# prettyPrintList
# properCase
# rotateVector
# scaleData
# scaleData
# solveODE
# variance
#
###################################################

###
# ACVStoPACS
##

"ACVStoPACS" <- function(acvs)
{
  checkVectorType(acvs,"numeric")
  p       <- length(acvs) - 1
  pacs    <- 1:p
  pacs[1] <- acvs[2]/acvs[1]

  if (p > 1) {

    phi <- pacs[1]
    var <- acvs[1] * (1 - pacs[1]^2)
    for(k in 2:p) {
      pacs[k] <- (acvs[k + 1] - sum(phi * rev(acvs[2:(length(phi) + 1)])))/var
      phi <- c(phi - pacs[k] * rev(phi), pacs[k])
      var <- var * (1 - pacs[k]^2)
    }
  }
  pacs
}

###
# aggregateData
###

"aggregateData" <- function(x, by, FUN, moving=FALSE, ...)
{
  x <- as.numeric(x)
  checkScalarType(by,"integer")
  checkRange(by, c(1,length(x)))
  if (!is.function(FUN))
    stop("FUN must be a function")

  n.sample <- length(x)

  if (!moving){
    n.usable <- (n.sample %/% by) * by

    z <- as.vector(apply(matrix(x[1:n.usable], nrow=by), MARGIN=2, FUN=FUN, ...))
    if (n.usable != n.sample)
      z <- c(z, as.vector(FUN(x[(n.usable+1):n.sample], ...)))
  }
  else{

  	moving  <- as.integer(moving)
  	checkRange(moving, c(1, n.sample))
    index   <- seq(moving)
    n.group <- as.integer(floor((n.sample - moving)/ by + 1))
    last    <- (n.group - 1)*by + moving

    if (!n.group)
      stop("Moving window not possible with input parameters")

    z <- vector("numeric", length=n.group)

    for (i in seq(n.group)){
      z[i] <- FUN(x[index + by*(i-1)], ...)
    }

    if (last < n.sample)
      z <- c(z, as.vector(FUN(x[(last+1):n.sample], ...)))
  }

  z
}

###
# checkRange
###

"checkRange" <- function(x, range.=0:1, inclusion=c(TRUE,TRUE)){
  checkVectorType(range.,"numeric")
  if (length(range.) != 2)
    stop("Range input must contain two elements")
  range. <- sort(range.)
  checkVectorType(inclusion,"logical")
  inclusion <- ifelse1(length(inclusion) == 1, rep(inclusion,2),
    inclusion[1:2])

  xrange <- range(x)
  left  <- ifelse1(inclusion[1], xrange[1] >= range.[1], xrange[1] > range.[1])
  right <- ifelse1(inclusion[2], xrange[2] <= range.[2], xrange[2] < range.[2])

  if (!(left && right))
    stop("Range of x not on specified interval")
  invisible(NULL)
}

###
# checkScalarType
###

"checkScalarType" <- function(x, isType="numeric"){

  if (!is.character(isType))
    stop("isType must be an object of class character")

  # In R, class(5) is "numeric", not "integer" so we
  # have to accommodate
  #
  # also, if S-PLUS in R parse mode, then integers are converted to floats (1 to 1.)
  # so we accommodate that as well. (set.parse.mode(-1) == 1) means that S-PLUS is in
  # "R" parse mode
  #if ((is.R() || (!is.R() && as.logical(set.parse.mode(-1)))) && isType == "integer"){
  if (isType == "integer")
  {
    if (!is.numeric(x) || length(x) > 1)
      stop(deparseText(substitute(x)), " must be scalar of class ", isType)
  }
  else{
    if (!eval(parse(text=paste("is.", isType, "(x)", sep=""))) || length(x) > 1)
      stop(deparseText(substitute(x)), " must be scalar of class ", isType)
  }

  invisible(NULL)
}

###
# checkVectorType
###

"checkVectorType" <- function(x, isType="numeric"){
  checkScalarType(isType,"character")

  # In R, class(c(1,3:5) is "numeric", not "integer" so we
  # have to accommodate
  #
  # also, if S-PLUS in R parse mode, then integers are converted to floats (1 to 1.)
  # so we accommodate that as well. (set.parse.mode(-1) == 1) means that S-PLUS is in
  # "R" parse mode
  #if ((is.R() || (!is.R() && as.logical(set.parse.mode(-1)))) && isType == "integer"){
  if (isType == "integer")
  {

    if (!isVectorAtomic(x) || !is.numeric(x))
      stop(deparseText(substitute(x)), " must be a vector of class ", isType)
  }
  else{

    if (!isVectorAtomic(x) || !eval(parse(text=paste("is.", isType, "(x)", sep=""))))
      stop(deparseText(substitute(x)), " must be a vector of class ", isType)
  }

  invisible(NULL)
}

###
# decibel
###

"decibel" <- function(x, type=1, na.zero=TRUE)
{
  checkScalarType(type,"integer")
  checkScalarType(na.zero,"logical")
  if (type < 1 || type > 2)
    stop("type must be on [1,2]")
  if (min(x) < 0)
   stop("cannot convert negative values to decibels")

  if (na.zero){
	izero <- which(x == 0)
	if (length(izero)){
      #warning("Zeros in data. Setting to NA before taking the logarithm")
      x[izero] <- NA
	}
  }

  ifelse1(type == 1, 10 * log10(x), 20 * log10(x))
}

###
# em
###

"em" <- function()
  c(strwidth("m"), strheight("m"))

###
# ilogb
###

"ilogb" <- function(x, base=2, eps=.Machine$double.eps * 1e9)
  as.integer(logb(x, base=base) + eps)

###
# isVectorAtomic
###

"isVectorAtomic" <- function(x)
  return(is.atomic(x) & any(c(NROW(x),NCOL(x)) == 1))

###
# linearSegmentation
###

"linearSegmentation" <- function(x, y, n.fit=5, angle.tolerance=5,
  aspect=TRUE, plot=FALSE, add=FALSE, ...)
{
  checkVectorType(x,"numeric")
  checkVectorType(y,"numeric")
  checkScalarType(n.fit,"integer")
  checkScalarType(angle.tolerance,"numeric")
  checkScalarType(aspect,"logical")

  if (angle.tolerance <= 0 || angle.tolerance >= 180)
    stop("angle.tolerance must be on the interval (0,180)")
  if (length(x) != length(y))
    stop("x and y must be the same length")
  if (n.fit < 3 || n.fit > length(x))
    stop("n.fit must be on the interval [3, length(x)]")

  # get rid of NAs
  bad <- union(which(is.na(y)), which(is.na(x)))

  if (length(bad) > 0){
    x <- x[-bad]
    y <- y[-bad]
  }

  # reverse if x not in increasing monotonic order
  if (any(diff(x) < 0))
  {
    x <- rev(x)
    y <- rev(y)
  }

  # initialize variables
  smooth   <- supsmu(x, y)$y
  y        <- smooth
  n.sample <- length(x)

  aspect.ratio    <- ifelse1(aspect, diff(range(x)) / diff(range(y)), 1)
  angle.tolerance <- angle.tolerance / aspect.ratio

  # initialize break vector
  breaks <- as.vector(itCall("RS_fractal_piecwise_linear_segmentation",
    matrix(as.double(x)), matrix(as.double(y)), n.fit, angle.tolerance)
    #COPY=rep(FALSE,4),
    #CLASSES=c("matrix", "matrix", "integer", "numeric"),
    #PACKAGE="ifultools")
    ) + 1

  attr(breaks, "smooth") <- y

  if (plot){
	  old.mfrow <- par(mfrow=c(1,1))
	  on.exit(par(old.mfrow))
	  plot(x=x, y=y, add=add)
  	abline(v=x[breaks], lty=2, ...)
  }

  breaks
}

###
# linearFit
###

"linearFit" <- function(x, y, fit=lmsreg, method="widest",
  n.fit=5, angle.tolerance=5, aspect=TRUE)
{
  checkVectorType(x,"numeric")
  checkVectorType(y,"numeric")
  if (!is.function(fit))
    stop("fit must be a linear regression function")
  checkScalarType(method,"character")
  checkScalarType(n.fit,"integer")
  checkScalarType(angle.tolerance,"numeric")
  checkScalarType(aspect,"logical")

  if (angle.tolerance <= 0 || angle.tolerance >= 180)
    stop("angle.tolerance must be on the interval (0,180)")
  if (length(x) != length(y))
    stop("x and y must be the same length")
  if (n.fit < 3 || n.fit > length(x))
    stop("n.fit must be on the interval [3, length(x)]")

  method <- match.arg(method, c("first", "last", "widest", "all"))

  # remove NAs
  nax <- which(is.na(x))
  nay <- which(is.na(y))
  bad <- union(nax, nay)

  if (length(bad) > 0){
    x <- x[-bad]
    y <- y[-bad]
  }

  N <- length(x)

  # obtain the appropriate scaling range
  # choose the range of scale which corresponds
  # to the smallest scales
  breaks <- linearSegmentation(x, y, n.fit=n.fit,
	angle.tolerance=angle.tolerance, aspect=aspect)

  nbreak <- length(breaks)

  if (!nbreak)
    method <- "all"

  if (method == "first")
    div <- seq(1, breaks[1])
  else if (method == "last")
    div <- seq(breaks[max(nbreak-1,1)], N)
  else if (method == "all")
    div <- seq(N)
  else if (method == "widest"){
    newbreaks <- unique(c(1,breaks,N))
    ib        <- order(diff(newbreaks))[length(newbreaks)-1]
    div       <- seq(newbreaks[ib], newbreaks[ib+1])
  }

  if (length(div) > 5)
    model <- fit(y ~ x, subset=div)
  else
    model <- lm(y ~ x, subset=div)

  return(model)
}

###
# logScale
###

"logScale" <- function(scale.min, scale.max,
  scale.ratio=2, scale.res=NULL, coerce=NULL)
{
  # check input arguments
  checkScalarType(scale.min, "numeric")
  checkScalarType(scale.max, "numeric")
  if (!is.null(scale.res))
    scale.ratio <- 1 / scale.res + 1
  checkScalarType(scale.ratio, "numeric")

  # check input arguments
  if (scale.ratio == 1)
    stop("scale.ratio cannot be unity")
  if (scale.ratio <= 0)
    stop("scale.ratio must be positive")

  n.scale <- ifelse1(scale.ratio > 1, ceiling(log(scale.max/scale.min)/log(scale.ratio)),
    ceiling(log(scale.min/scale.max)/log(scale.ratio)))

  scale    <- vector("numeric", length=n.scale)
  scale[1] <- ifelse1(scale.ratio > 1, scale.min, scale.max)

  for (i in seq(2, n.scale)){
    scale[i] <- scale[i - 1] * scale.ratio
  }

  if (!is.null(coerce) && is.function(coerce))
    scale <- coerce(scale)

  unique(scale)
}

###
# mergeList
###

"mergeList" <- function(x,y)
{
  # given lists x and y, replace
  # the common named objects in x with
  # those of y and append the uncommon y components
  # to x and return

  if (!is.list(y))
    stop("y must be a list")
  if (is.null(x))
    return(y)
  if (!is.list(x))
    stop("x must be a list")
  x[names(y)] <- y
  x
}

###
# nextDyadic
###

"nextDyadic" <- function(x)
{
  # returns the smallest integer n such that
  # (i)  n is a power of two (e.g., 1, 2, 4, 8, 16, ...) and
  # (ii) n is greater than or equal to x
  x[x < 1] <- 1
  as.integer(round(2^(ceiling(logb(x, base=2)))))
}

###
# ordinal
###

"ordinal" <- function(x){
  checkScalarType(x,"numeric")
  z <- strsplit(as.character(round(x)),"")[[1]]

  paste(paste(z, collapse=""), switch(as.integer(z[length(z)])+1,
    "th","st","nd","rd","th","th","th","th","th","th"),sep="")
}

###
# prettyPrintList
###

"prettyPrintList" <- function(x, header=NULL, justify="left", sep=":")
{
  if (!is.list(x) || is.null(names(x)))
    stop("x must be a list containing named objects")
  if (!is.null(header) && (!is.character(header) || length(header) > 1))
    stop("header must be a single character string")
  if (!is.character(justify) || length(justify) > 1)
    stop("justify must be a single character string")
  if (!is.character(sep) || length(sep) > 1)
    stop("sep must be a single character string")

  justify <- match.arg(justify, c("none","left","right","decimal"))

  if (!is.null(header))
    cat(header,"\n", rep("-",nchar(header)),"\n",sep="")

  # prune list of NULL values.
  # if x <- list("some really large name"=NULL, cat="dog")
  # the spearator will be spaced way too far to the right
  # due to the influence of the first entry name. thus,
  # we eliminate any such NULL entries altogether to
  # avoid this problem
  x <- x[!unlist(lapply(x, is.null))]

  if (!length(x))
    return(invisible(NULL))

  categories <- format(names(x), justify=justify)

  for (i in seq(along=categories)){
    if (!is.null(x[[i]]))
      cat(categories[i], sep, x[[i]], "\n", sep=" ")
  }

  invisible(NULL)
}

###
# properCase
###

"properCase" <- function(x, pre="", post="", sep="")
{
  if (!is.character(x))
    stop("Input must be of class character")
  checkScalarType(pre, "character")
  checkScalarType(post,"character")
  checkScalarType(sep, "character")

  # make first letter uppercase
  substring(x,1,1) <- upperCase(substring(x,1,1))

  # make remainder lower case
  substring(x,2) <- lowerCase(substring(x,2))

  # add pre and post strings
  z <- paste(pre,x,post,sep=sep)

  # maintain name if available
  names(z) <- names(x)

  z
}

###
# rotateVector
###

"rotateVector" <- function(x, shift=1)
{
  # Rotate (circularly wrap) a vector by "shift" units to the right.
  # If shift < 0, the data is permuted to the left by abs(shift) units.
  N <- length(x)

  # reduce shifts beyond period one
  shift <- sign(shift) * (abs(shift) %% N)

  if (shift < 0)
    x[c((abs(shift) + 1):N, 1:abs(shift))]
  else if (shift > 0)
    x[c((N - shift + 1):N, 1:(N - shift))]
  else
    x
}

###
# scaleData
###

"scaleData" <- function(x, scale.="linear")
{
  checkScalarType(scale.,"character")
  supported.scales <- c("linear","log2","log10","log","db")
  scale. <- match.arg(lowerCase(scale.), supported.scales)

  z <- switch(scale.,
    linear=x,
    log2  =logb(x, base=2),
    log10 =logb(x, base=10),
    log   =log(x),
    db    =decibel(x))

  attributes(z) <- c(attributes(x), list(scalestr=switch(scale.,linear="",log2="log2",log10="log10",log="log",db="dB")))
  z
}

###
# solveODE
###

"solveODE" <- function(FUN, initial=NULL, step=0.01, n.sample=1000, n.transient=100, ...){

  if (is.null(initial))
    initial <- FUN(x=NULL)

  if (n.transient < 0)
    stop("Number of transient points must be positive or zero")

  x <- initial

  result <- matrix(NA, nrow=n.sample + n.transient, ncol=length(x))

  result[ 1, ] <- x

  for (n in seq(2, n.sample + n.transient)){

    ## calculate four approximations
    m1 <- step * FUN(x, ...)
    m2 <- step * FUN(x + m1 / 2, ...)
    m3 <- step * FUN(x + m2 / 2, ...)
    m4 <- step * FUN(x + m3, ...)

    # estimate the next iterate via (fourth order)
    # Runge-Kutta approximation
    xnew <- x + (m1 + 2 * m2 + 2 * m3 + m4) / 6

    # record the result
    result[ n, ] <- xnew

    # update the current state
    x <- xnew
  }

  # remove the transient and return result
  return(data.frame(result[-seq(n.transient), ]))
}

###
# variance
###

"variance" <- function(x, na.rm=TRUE, unbiased=FALSE)
{
  if (is(x,"signalSeries"))
    x <- as(x,"numeric")

  N <- numRows(x)

  v <- var(x,na.rm=na.rm)
  if (!unbiased) v <- v*(N-1)/N

  v
}
