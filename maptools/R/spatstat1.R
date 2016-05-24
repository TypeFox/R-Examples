# as.ppp method to be used in spatstat:

as.ppp.SpatialPoints = function(X) {
	# require(spatstat)
	if (!requireNamespace("spatstat", quietly = TRUE))
		stop("package spatstat required for as.ppp.SpatialPoints")
        bb <- bbox(X)
        colnames(bb) <- NULL
	W = spatstat::owin(bb[1,], bb[2,])
	cc = coordinates(X)
	spatstat::ppp(cc[,1], cc[,2], window = W, marks = NULL, check=FALSE)
}

setAs("SpatialPoints", "ppp", function(from) as.ppp.SpatialPoints(from))

# Mike Sumner 20101011
as.ppp.SpatialPointsDataFrame = function(X) {
	# require(spatstat)
	if (!requireNamespace("spatstat", quietly = TRUE))
		stop("package spatstat required for as.ppp.SpatialPointsDataFrame")
	bb <- bbox(X)
        colnames(bb) <- NULL
	W  <- spatstat::owin(bb[1,], bb[2,])
	nc <- ncol(X)
        marks <- if(nc == 0) NULL else slot(X, "data")
#        if(nc > 1)
#          warning(paste(nc-1, "columns of data frame discarded"))
	cc <- coordinates(X)
	spatstat::ppp(cc[,1], cc[,2], window = W, marks = marks, check=FALSE)
}

setAs("SpatialPointsDataFrame", "ppp", function(from) as.ppp.SpatialPointsDataFrame(from))

as.owin.SpatialGridDataFrame = function(W, ..., fatal) {
	# require(spatstat)
	if (!requireNamespace("spatstat", quietly = TRUE))
		stop("package spatstat required for as.owin.SpatialGridDataFrame")
	# W = from
	m = t(!is.na(as(W, "matrix")))
	spatstat::owin(bbox(W)[1,], bbox(W)[2,], mask = m[nrow(m):1,])
}

setAs("SpatialGridDataFrame", "owin", function(from) as.owin.SpatialGridDataFrame(from))

as.owin.SpatialPixelsDataFrame = function(W, ..., fatal) {
	# require(spatstat)
	if (!requireNamespace("spatstat", quietly = TRUE))
		stop("package spatstat required for as.owin.SpatialPixelsDataFrame")
	# W = from
	m = t(!is.na(as(W, "matrix")))
	spatstat::owin(bbox(W)[1,], bbox(W)[2,], mask = m[nrow(m):1,])
}

setAs("SpatialPixelsDataFrame", "owin", function(from) as.owin.SpatialPixelsDataFrame(from))

as.owin.SpatialPolygons = function(W, ..., fatal) {
	# require(spatstat)
	# W = from
	if (!inherits(W, "SpatialPolygons")) 
		stop("W must be a SpatialPolygons object")
	res <- .SP2owin(W)
	res
}

setAs("SpatialPolygons", "owin", function(from) as.owin.SpatialPolygons(from))

# methods for coercion to Spatial Polygons by Adrian Baddeley

owin2Polygons <- function(x, id="1") {
	if (!requireNamespace("spatstat", quietly = TRUE))
		stop("package spatstat required for as.owin.SpatialPixelsDataFrame")
  stopifnot(spatstat::is.owin(x))
  x <- spatstat::as.polygonal(x)
  closering <- function(df) { df[c(seq(nrow(df)), 1), ] }
  if (x$type == "polygonal") {
      pieces <- lapply(x$bdry,
                   function(p) {
                     Polygon(coords=closering(cbind(p$x,p$y)),
                             hole=spatstat::is.hole.xypolygon(p))  })
  } else if (x$type == "rectangle") {
      rectCrds <- cbind(x$xrange[c(1,1,2,2,1)], x$yrange[c(1,2,2,1,1)])
      pieces <- list(Polygon(rectCrds, hole=FALSE))
  } else stop("owin2Polygons: unknown type:", x$type)
  z <- Polygons(pieces, id)
  return(z)
}

as.SpatialPolygons.tess <- function(x) {
	#require(spatstat)
	if (!requireNamespace("spatstat", quietly = TRUE))
		stop("package spatstat required for as.owin.SpatialPixelsDataFrame")
  stopifnot(spatstat::is.tess(x))
  y <- spatstat::tiles(x)
  nam <- names(y)
  z <- list()
  for(i in seq(y)) {
    zi <- try(owin2Polygons(y[[i]], nam[i]), silent=TRUE)
    if (class(zi) == "try-error") {
      warning(paste("tile", i, "defective\n", as.character(zi)))
    } else {
      z[[i]] <- zi
    }
  }
  return(SpatialPolygons(z))
}

setAs("tess", "SpatialPolygons", function(from) as.SpatialPolygons.tess(from))


as.SpatialPolygons.owin <- function(x) {
	if (!requireNamespace("spatstat", quietly = TRUE))
		stop("package spatstat required for as.owin.SpatialPixelsDataFrame")
	#require(spatstat)
  stopifnot(spatstat::is.owin(x))
  y <- owin2Polygons(x)
  z <- SpatialPolygons(list(y))
  return(z)
}

setAs("owin", "SpatialPolygons", function(from) as.SpatialPolygons.owin(from))



# methods for 'as.psp' for sp classes by Adrian Baddeley

as.psp.Line <- function(from, ..., window=NULL, marks=NULL, fatal) {
	#require(spatstat)
	if (!requireNamespace("spatstat", quietly = TRUE))
		stop("package spatstat required for as.owin.SpatialPixelsDataFrame")
  xy <- slot(from, "coords")
  df <- as.data.frame(cbind(xy[-nrow(xy), , drop=FALSE], xy[-1, ,
drop=FALSE]))
  if(is.null(window)) {
    xrange <- range(xy[,1])
    yrange <- range(xy[,2])
    window <- spatstat::owin(xrange, yrange)
  }
  return(spatstat::as.psp(df, window=window, marks=marks))
}

setAs("Line", "psp", function(from) as.psp.Line(from))
  
as.psp.Lines <- function(from, ..., window=NULL, marks=NULL, fatal) {
	# require(spatstat)
	if (!requireNamespace("spatstat", quietly = TRUE))
		stop("package spatstat required for as.owin.SpatialPixelsDataFrame")
  y <- lapply(slot(from, "Lines"), as.psp.Line, window=window)
  if (as.character(packageVersion("spatstat")) < "1.22.0") {
    z <- spatstat::superimposePSP(y, window=window)
  } else {
    z <- do.call(spatstat::superimpose,c(y,list(W=window)))
#    z <- superimpose(y, window=window)
  }
  if(!is.null(marks))
    spatstat::marks(z) <- marks
  return(z)
}

setAs("Lines", "psp", function(from) as.psp.Lines(from))

as.psp.SpatialLines <- function(from, ..., window=NULL, marks=NULL,
                                 characterMarks=FALSE, fatal) {
	# require(spatstat)
	if (!requireNamespace("spatstat", quietly = TRUE))
		stop("package spatstat required for as.owin.SpatialPixelsDataFrame")
  if(is.null(window)) {
    w <- slot(from, "bbox")
    window <- spatstat::owin(w[1,], w[2,])
  }
  lin <- slot(from, "lines")
  y <- lapply(lin, as.psp.Lines, window=window)
  id <- row.names(from)
  if(is.null(marks))
    for (i in seq(y)) 
      spatstat::marks(y[[i]]) <- if(characterMarks) id[i] else factor(id[i])
# modified 110401 Rolf Turner
#    for(i in seq(y)) 
#      marks(y[[i]]) <- id[i]
  if (as.character(packageVersion("spatstat")) < "1.22.0") {
    z <- do.call(spatstat::superimposePSP, list(y, window=window))
  } else {
#    z <- do.call("superimpose", list(y, window=window))
    z <- do.call(spatstat::superimpose, c(y, list(W = window)))
  }
  if(!is.null(marks))
    spatstat::marks(z) <- marks
  return(z)
}

setAs("SpatialLines", "psp", function(from) as.psp.SpatialLines(from))

as.psp.SpatialLinesDataFrame <- function(from, ..., window=NULL, marks=NULL, fatal) {
	# require(spatstat)
	if (!requireNamespace("spatstat", quietly = TRUE))
		stop("package spatstat required for as.owin.SpatialPixelsDataFrame")
  y <- as(from, "SpatialLines")
  z <- spatstat::as.psp(y, window=window, marks=marks)
  if(is.null(marks)) {
    # extract marks from first column of data frame
    df <- from@data
    if(is.null(df) || (nc <- ncol(df)) == 0)
      return(z)
    if(nc > 1) 
      warning(paste(nc-1, "columns of data frame discarded"))
    marx <- df[,1]
    nseg.Line  <- function(x) { return(nrow(x@coords)-1) }
    nseg.Lines <- function(x) { return(sum(unlist(lapply(x@Lines, nseg.Line)))) }
    nrep <- unlist(lapply(y@lines, nseg.Lines))
    spatstat::marks(z) <- rep(marx, nrep)
  }
  return(z)
}

setAs("SpatialLinesDataFrame", "psp", function(from) as.psp.SpatialLinesDataFrame(from))

# 111117 from psp to SpatialLines, Rolf Turner, Adrian Baddeley, Mathieu Rajerison

as.SpatialLines.psp <- function(from) {

     ends2line <- function(x) Line(matrix(x, ncol=2, byrow=TRUE))
     munch <- function(z) { Lines(ends2line(as.numeric(z[1:4])), ID=z[5]) }
  
     ends <- as.data.frame(from)[,1:4]
     ends[,5] <- row.names(ends)
     y <- apply(ends, 1, munch)
     SpatialLines(y)
}

setAs("psp", "SpatialLines", function(from) as.SpatialLines.psp(from))
