################################################
## splus2R series functionality
##
##  Constructor Functions and methods:
##
##    numericSequence
##
##    seriesData
##
##      asSeriesData
##      seriesData<-
##      seriesDataNew
##      seriesDataValid
##
##    signalSeries
##
##      as.data.frame.signalSeries
##      as.matrix.signalSeries
##      cumsum.signalSeries
##      deltat.signalSeries
##      diff.signalSeries
##      plot.signalSeries
##
################################################

###
# numericSequence
###

"numericSequence" <- function(from, to, by, length.)
{
  # creation function for sequences
  # no args -> create default object
  if(missing(from) && missing(to) && missing(by) && missing(length.))
    return(new("numericSequence"))
  if(missing(from))
    from <- numeric(0)
  if(missing(to))
    to <- numeric(0)
  if(missing(by))
    by <- numeric(0)
  if(missing(length.))
    length. <- numeric(0)
  lens <- c(length(from), length(to), length(by), length(length.))
  if(any(lens > 1) || (sum(lens) != 3))
    stop(paste(
      "Exactly three of from, to, by, length. must be specified",
      "and they must each be of length 1"))
  if(anyMissing(c(from, to, by, length.)))
    stop("NA values not allowed in sequences")
  ret <- new("numericSequence", from = as(from, "numeric"), to = as(
    to, "numeric"), by = as(by, "numeric"), length = as(length.,
    "integer"))
  ret
}

################################################
##
##    seriesData
##
##      asSeriesData
##      seriesData<-
##      seriesDataNew
##      seriesDataValid
##
################################################

###
# seriesData
###

"seriesData" <- function(object)
{
	# return the data inside an ordered data object
	object@data
}

###
# asSeriesData
###

"asSeriesData" <- function(x)
  as.rectangular(x)

###
# seriesData<-
###

"seriesData<-" <- function(object, value)
{
  # replace the data inside an ordered data object
  value <- asSeriesData(value)
  if (length(object@positions) != numRows(value))
    stop("Positions and data lengths do not agree")
  object@data <- value
  object
}

###
# seriesDataNew
###

"seriesDataNew" <- function(){
  # R does not support matrix(NULL), i.e.,
  # NULL matrices are apparently not allowed
  # just use matrix(numeric()) instead

  #  matrix(NULL)
  matrix(numeric())
}

###
# seriesDataValid
###

"seriesDataValid" <- function(object)
  is.rectangular(object)

################################################
##
##    signalSeries
##
##      as.data.frame.signalSeries
##      as.matrix.signalSeries
##      cumsum.signalSeries
##      deltat.signalSeries
##      diff.signalSeries
##      plot.signalSeries
##
################################################

###
# signalSeries
###

"signalSeries" <- function(data, positions., units, units.position, from=1, by=1)
{
  # function to create a signalSeries object
  # positions, if supplied, overrides from, to, by
  if (missing(positions.) && missing(data) && missing(units) &&
    missing(units.position) && missing(from) && missing(by))
      return(new("signalSeries"))

  if (!missing(positions.) && (length(positions.) != numRows(data)))
    stop("Positions and data lengths do not agree")

  ret <- new("signalSeries")
  ret@data <- asSeriesData(data)

  if (missing(positions.)) {
    len <- numRows(ret@data)
    ret@positions <- numericSequence(from = from, length. = len, by = by)
  }
  else{
    if (!is(positions., "numericSequence"))
      positions. <- as(positions., "numericSequence")
    ret@positions <- positions.
  }

  if (!missing(units))
    ret@units <- as(units, "character")

  if (!missing(units.position))
    ret@units.position <- as(units.position, "character")
  ret
}

###
# as.data.frame.signalSeries
###

"as.data.frame.signalSeries" <-
  function (x, row.names = NULL, optional = FALSE, ...) 
{
  as.data.frame(x@data,
                row.names = if(is.null(row.names))
                as(x@positions, "character") else row.names,
                optional = optional, ...)
}

###
# as.matrix.signalSeries
###

"as.matrix.signalSeries" <- function(x, ...)
  as.matrix(x@data, ...)

###
# cumsum.signalSeries
###

"cumsum.signalSeries" <- function(x)
{
   data <- cumsum(x@data)
   x@data <- data
   x
}

###
# deltat.signalSeries
###

"deltat.signalSeries" <- function(x, ...) positions(x)@by

###
# diff.signalSeries
###

"diff.signalSeries" <- function(x,...)
{
   data <- diff(x@data,...)
   pos  <- as(positions(x)[-1],"numericSequence")
   signalSeries(data, positions.=pos, units=x@units, units.position=x@units.position)
}

###
# plot.signalSeries
###

"plot.signalSeries" <-
  function(x, y, ..., main=NULL, ylab=x@units[1], xlab=x@units.position,
  top.ticks = FALSE, right.ticks = FALSE, reference.grid=TRUE, merge.args = list(
  pos = "union", how = "interp"), x.axis.args = list(), y.axis.args =
  list(), plot.args = list(), log.axes = "", complex.convert = Mod,
  dB = FALSE, frame = sys.nframe(), col=NULL, type="l")
{
  mergeArgList <- function (x, y) {
    if (!is.list(y))
        stop("y must be a list")
    if (is.null(x))
        return(y)
    if (!is.list(x))
        stop("x must be a list")
    x[names(y)] <- y
    x
  }

  y    <- c(list(x), if(!missing(y)) list(y), list(...))
  ny   <- length(y)

  ylim <- range(unlist(lapply(y, function(x) range(seriesData(x)))))
  xlim <- range(unlist(lapply(y, function(x) range(as(positions(x),"numeric")))))
  plot(xlim, ylim, type="n", xlab=xlab, ylab=ylab, log=log.axes)

  if (is.null(col))
    col <- seq(ny)
  else
    col <- rep(col, length = ny)

  type <- rep(type, length = ny)

  if (is.null(main))
    main <- paste(unlist(lapply(y,function(x) x@title)), sep=",")

  for (i in seq(along=y)){
    xdata <- as(positions(y[[i]]),"numeric")
    ydata <- y[[i]]@data
    do.call("lines",
	    mergeArgList(list(x=xdata, y=ydata, col=col[i], type=type[i]),
			 plot.args))
    # mergeArgList is needed in R, not S+
    # (S+ allows duplicate args, uses last)
  }

  title(main=main)
  invisible(NULL)
}










