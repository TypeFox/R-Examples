###
# Make functions S4 generic
###

#setGeneric("as.data.frame")
setGeneric("as.numeric")
setGeneric("as.vector")
setGeneric("duplicated")
setGeneric("match")
setGeneric("mean")
setGeneric("median")
setGeneric("ncol")
#setGeneric("nrow")
setGeneric("plot")
setGeneric("quantile")
setGeneric("rev")
setGeneric("sort")
#setGeneric("sort")
#setGeneric("sort.list")
#setGeneric("sort.list")
setGeneric("unique")
setGeneric("which.na")

###
# Class: numericSequence
###

setClass("numericSequence",
  representation(
    from = "numeric",
    to = "numeric",
    by = "numeric",
    length = "integer"),
  prototype=prototype(from = 0, to = 0, by = numeric(0), length = as.integer(0)))

setAs("numericSequence", "numeric",
  function(from){
    if(!length(from))
      return(numeric(0))

    # must have exactly one of the slot vectors missing, or
    # else we will ignore the length
    # call seq appropriately

    if(is.missing(from@from))
      seq(to = from@to, by = from@by, length = from@length)
    else if(is.missing(from@to))
      seq(from = from@from, by = from@by, length = from@length)
    else if(is.missing(from@by))
      seq(from = from@from, to = from@to, length = from@length)
    else{
      if(!is.missing(from@length))
        warning("Ignoring length slot in sequence")
      seq(from = from@from, to = from@to, by = from@by)
    }
  })

setAs("numericSequence", "integer",
  function(from) as(as(from, "numeric"), "integer"))

setAs("numericSequence", "character",
  function(from) as(as(from, "numeric"), "character"))

setAs("numeric", "numericSequence",
  function(from){

    # convert numeric vector to a sequence, if possible
    len <- length(from)
    if(!len)
      return(numericSequence())
    if(len == 1)
      return(numericSequence(from, from, length. = 1))

    # verify that it's a valid sequence
    diffs <- from[ -1 ] - from[ - length(from) ]
    tol   <- options("sequence.tol")[[1]]
    if(any(abs(diffs - diffs[1]) > tol))
      stop("Numeric vector is not a regularly spaced sequence")

    # create sequence object
    numericSequence(from = from[1], by = diffs[1],
      length. = len)
  })

setMethod("[", signature(x = "numericSequence"),
  function(x, i, j, ..., drop = TRUE){
    x <- as(x, "numeric")
    if (missing(i))
      i <- seq(x)
    callGeneric(x, i, ..., drop = TRUE)
  })

setMethod("[[", signature(x = "numericSequence"),
  function(x, i, j, ...) {
    x <- as(x, "numeric")
    if (missing(i))
      i <- seq(x)
    callGeneric(x, i, ..., drop = TRUE)
  })

setReplaceMethod("[", signature(x = "numericSequence"),
  function(x, i, j, ..., value) {
     x <- as(x, "numeric")
     x[...] <- value
     x
  })

setReplaceMethod("[[", signature(x = "numericSequence"),
  function(x, i, j, ..., value){
    x <- as(x, "numeric")
    x[[...]] <- value
    x
   })

setMethod("length", signature(x = "numericSequence"),
  function(x){

    # length for sequences
    slots.missing <- c(is.missing(x@from),
      is.missing(x@to),
      is.missing(x@by),
      is.missing(x@length))
      how.many.there <- sum(c(1,1,1,1)[!slots.missing])

    # if it has a length slot, and missing one of others, that's the length
    if(!slots.missing[4] && (how.many.there < 4))
       return(x@length)

    # if the by slot is 0, then from must be
    # equal to to and it's a length 1 sequence
    # assuming it's valid!
    if(!x@by)
      return(1)
    # calculate the length from to/from/by
    as.integer(floor(1 + ((x@to - x@from) / x@by)))
  })

setMethod("show", "numericSequence", function(object)
{
  if(!is.missing(object@from))
    cat("from:  ", object@from, "\n")
  if(!is.missing(object@to))
    cat("to:    ", object@to, "\n")
  if(!is.missing(object@by))
    cat("by:    ", object@by, "\n")
  if(!is.missing(object@length))
    cat("length:", object@length, "\n")
  nums <- format(as(object,"numeric"))
  if(length(nums) > 5)
    nums <- c(nums[1:3], "...", nums[length(nums)])
  invisible(print(nums, quote=FALSE))
})

setMethod("summary", "numericSequence", function(object, ...)
{
  sumry <- numeric(0)
  nms <- character(0)

  from <- object@from
  to <- object@to
  by <- object@by
  len <- object@length

  if(is.missing(from))
    from <- to - (len - 1) * by
  if(is.missing(to))
    to <- from + (len - 1) * by
  if(is.missing(by))
    by <- (to - from) / (len - 1)
  if(is.missing(len))
    len <- floor((to - from) / by) + 1

  sumry <- matrix(c(from, to, by, len), nrow = 1,
       dimnames = list("", c("From", "To", "By", "Length")))
  oldClass(sumry) <- "table"
  sumry
})

setMethod("Math", "numericSequence",
    function(x) callGeneric(as(x, "numeric")))

setMethod("Math2", signature(x = "numericSequence"),
    function(x, digits) callGeneric(as(x, "numeric"), digits))

setMethod("Summary", signature(x = "numericSequence"),
    function(x, ..., na.rm= FALSE)
            callGeneric(as(x, "numeric"), ..., na.rm = na.rm))

setMethod("Ops", signature(e1 = "numericSequence"),
    function(e1, e2 = NULL)
      callGeneric(as(e1, "numeric"), e2))
setMethod("Ops", signature(e2 = "numericSequence"),
    function(e1, e2 = NULL)
      callGeneric(e1, as(e2, "numeric")))

setMethod("which.na", "numericSequence",
  function(x){
    # valid sequence objects cannot have NA
    numeric(0)
  })

setMethod("is.na", "numericSequence",
  function(x){
    # valid sequence objects cannot have NA
    rep(FALSE,length(x))
  })

setMethod("match", signature(x = "numericSequence"),
  function(x, table, nomatch = NA, incomparables= FALSE)
    match(as(x, "numeric"), table, nomatch, incomparables))

setMethod("match", signature(table = "numericSequence"),
  function(x, table, nomatch = NA, incomparables= FALSE)
    match(x, as(table, "numeric"), nomatch, incomparables))

setMethod("unique", signature(x = "numericSequence"),
           function(x, incomparables = FALSE, ...){
	     ifelse1(identical(x@from, x@to),
		     x@from,
		     length(x@by) && x@by == 0,
		     c(x@from, x@to), # one of from & to is numeric(0)
		     x)
	   })

setMethod("duplicated", signature(x = "numericSequence"),
           function(x, incomparables= FALSE ){
	     ifelse1(identical(x@from, x@to) ||
		     (length(x@by) && x@by == 0),
		     ifelse1(x@length,
			     c(FALSE, rep(TRUE, x@length-1)),
			     rep(FALSE, x@length)),
		     rep(FALSE, length(x)))
	   })

setMethod("mean", signature(x = "numericSequence"),
           function(x, trim = 0, na.rm= FALSE)
               mean(as(x, "numeric"), trim, na.rm))
setMethod("median", signature(x = "numericSequence"),
           function(x, na.rm= FALSE)
               median(as(x, "numeric"), na.rm))
setMethod("quantile", signature(x = "numericSequence"),
           function(x, probs = seq(0, 1, 0.25), na.rm= FALSE)
               quantile(as(x, "numeric"), probs, na.rm))

setMethod("rev", signature(x = "numericSequence"),
  function(x) {
    # Reverse a numericSequence object the cheap way - switch slots
    tmp <- x@from
    x@from <- x@to
    x@to <- tmp
    if(!is.missing(x@by))
      x@by <- -x@by
    x
  })


setMethod("sort", signature(x = "numericSequence"),
  function(x, decreasing = FALSE, ...){
    # sort a numericSequence object (reverse if necessary)
    if(is.missing(x@by)){
      if((!decreasing && x@to < x@from) ||
	 (decreasing && x@to > x@from))
        x <- rev(x)
    }
    else{
      if((!decreasing && x@by < 0) ||
	 (decreasing && x@by > 0))
	x <- rev(x)
    }
  })


#setMethod("sort.list", signature(x = "numericSequence"),
#  # I tried letting the argument list be (x, ...).  That failed.
#  function(x, partial = NULL, na.last = TRUE, decreasing = FALSE,
#    method = c("shell", "quick", "radix")){
#    # sort.list method for a numericSequence object
#    # ... may include decreasing - pass that to sort
#    sort(x, decreasing = decreasing)
#  })


###
# Class: signalSeries
###

setClass("signalSeries",
  representation(
    data           = "ANY",
    positions      = "numericSequence",
    units          = "character",
    title          = "character",
    documentation  = "character",
    attributes     = "ANY",
    units.position = "character"),
    prototype = prototype(
       data           = seriesDataNew(),
       positions      = numericSequence(),
       units          = character(),
       title          = character(),
       documentation  = character(),
       units.position = character())
)

setAs("signalSeries", "numeric",
      function(from) as(from@data, "numeric"))
setAs("signalSeries", "complex",
      function(from) as(from@data, "complex"))
setAs("signalSeries", "character",
      function(from) as(from@data, "character"))
setAs("signalSeries", "matrix",
      function(from) as(from@data, "matrix"))
setAs("signalSeries", "logical",
      function(from) as(from@data, "logical"))
setAs("signalSeries", "integer",
      function(from) as(from@data, "integer"))
setAs("signalSeries", "vector",
      function(from) as(from@data, "vector"))
setAs("signalSeries", "data.frame",
      function(from) as(from@data, "data.frame"))

setMethod("show", "signalSeries", function(object)
{
  # print out the data by making a table
  if(!numRows(object@data))
  {
    cat("NULL data\n")
    return()
  }

  if(is.data.frame(object@data))
    cols <- unlist(lapply(1:numCols(object@data),
       function(i, x)
         format(x[,i,drop=TRUE]),
         object@data))
  else
    cols <- format(as(object@data, "vector"))

  cols <- c(format(as(object@positions, "numeric")), cols)

  cols <- matrix(cols, ncol = numCols(object@data) + 1)
  colnames <- colIds(object@data)
  if(is.null(colnames) || !length(colnames))
    colnames <- 1:numCols(object@data)

  dimnames(cols) <- list(rep("",numRows(cols)),
         c("Positions", format(colnames)))

  oldClass(cols) <- "table"
  show(cols)
  NULL
})

setMethod("summary", "signalSeries", function(object, ...)
{
  ps <- summary(object@positions)
  # assume summary() always returns a 1-row table; want matrix
  oldClass(ps) <- NULL
  ds <- lapply(1:numCols(object@data),
          function(i, x)
          {
      ret <- summary(subscript2d(x,,i))
            oldClass(ret) <- NULL
            ret
    }, object@data)

  lens <- sapply(ds, "numCols")
  lens <- c(numCols(ps), lens)
  maxlen <- max(lens)

  # make all summaries into name : value character strings

  pastenames <- function(tab, newlen)
  {
    ret <- paste(format(dimnames(tab)[[2]]), ":",  format(tab),
      "  ", sep = "")
    length(ret) <- newlen
    ret
  }

  allcols <- cbind(pastenames(ps, maxlen),
        sapply(ds, pastenames, maxlen))

  colnames <- colIds(object@data)
  if(is.null(colnames))
    colnames <- 1:numCols(object@data)
  dimnames(allcols) <- list(rep("",maxlen),
         c("Positions", colnames))
  oldClass(allcols) <- "table"
  allcols
})

#setMethod("nrow", "signalSeries", function(x, columns) numRows(x@data))
#setMethod("numRows", "signalSeries", function(x, columns) numRows(x@data))
setMethod("ncol", "signalSeries", function(x) numCols(x@data))
#setMethod("numCols", "signalSeries", function(x) numCols(x@data))
setMethod("dim", "signalSeries", function(x) dim(x@data))
setMethod("length", "signalSeries", function(x) length(x@data))
#
#setGeneric("seriesLength", function(x) standardGeneric("seriesLength"))
#setMethod("seriesLength", "signalSeries", function(x) length(x@positions))
#setMethod("format", "signalSeries", function(x, ...) format(x@data, ...))

setMethod("[", "signalSeries", function(x, i, j, ..., drop = TRUE)
{
  if (missing(i))
    i <- seq(length(x))
  if (missing(j))
    j <- 1

  if (j > 1)
    stop("index j is out of range")
  if (i < 0 || i > length(x))
    stop("index i is out of range")

  pos  <- as(positions(x),"numeric")[i]
  data <- x@data[i]

  z <- signalSeries(data=data, positions.=as(pos,"numericSequence"),
    units.position=x@units.position)

  z@title <- x@title
  z@documentation <- x@documentation
  z
})

#setMethod("rowIds", "seriesVirtual", function(x) x@positions)
#setMethod("colIds", "seriesVirtual", function(x) colIds(x@data))

setMethod("plot", "signalSeries",
  function(x, y, ...) invisible(plot.signalSeries(x, y, ...)))

setMethod("as.vector", "signalSeries", function(x, mode="any") x@data)
setMethod("as.numeric", "signalSeries", function(x, mode="numeric") x@data)

#setMethod("as.data.frame","signalSeries",
#  function(x, mode="data.frame", row.names, optional)
#    data.frame(x@data))

setMethod("mean", signature(x ="signalSeries"),
   function(x, trim=0, na.rm=FALSE)
      mean(as(x, "numeric"), trim, na.rm))

setAs("list", "signalSeries",
  function(from){

    nms <- c("data","positions","start.position","end.position","future.positions","units",
       "title","documentation","units.position")

    if (!all(is.element(nms, names(from))))
      stop("Input list must contain the following named objects:\n  ",
        paste(nms, collapse=", "))

    pos <- from$positions
    pos <- numericSequence(from = pos$from, by = pos$by, to = pos$to, length. = pos$length)

    z <- signalSeries(from$data, positions.=pos, from$units, from$units.position)
    z@title <- from$title
    z@documentation <- from$documentation
    z
  }
)


setMethod("as.numeric", "signalSeries", function(x, mode="numeric") x@data)

###
# signalSeries Group Generic Functions
###

setMethod("Arith", "signalSeries", function(e1, e2){
  e1@data <- asSeriesData(callGeneric(e1@data, e2@data))
  e1
})

setMethod("Compare","signalSeries", function(e1, e2)
{
  e1@data <- asSeriesData(callGeneric(e1@data, e2@data))
  e1
})

setMethod("Ops", "signalSeries", function(e1, e2=NULL){
  e1@data <- asSeriesData(callGeneric(e1@data, e2))
  e1
})

setMethod("Logic","signalSeries", function(e1, e2)
{
  e1@data <- asSeriesData(callGeneric(e1@data, e2@data))
  e1
})

setMethod("Math","signalSeries",function(x)
{
  x@data <- asSeriesData(callGeneric(x@data))
  x
})

setMethod("Math2","signalSeries",function(x, digits)
{
  x@data <- asSeriesData(callGeneric(x@data, digits))
})

## There is no dotInternal(sum()) in R: it is primitive!
## Altered in an NMU 2012-09-21

## TODO: The following doesn't work for some reason
setMethod("Summary","signalSeries", function(x, ..., na.rm=FALSE)
  sum(..., na.rm = na.rm))

setMethod("min","signalSeries", function(x, ..., na.rm=FALSE)
  min(as(x,"numeric"), ..., na.rm=na.rm))

setMethod("sum","signalSeries", function(x, ..., na.rm=FALSE)
  sum(x@data, ..., na.rm = na.rm))

