#' Function to handle animal track data, organized as \code{"trip"}s
#'
#'
#' Create an object of class \code{"trip"}, extending the basic functionality
#' of \code{\link[sp]{SpatialPointsDataFrame}} by specifying the data columns
#' that define the "TimeOrdered" quality of the records.
#'
#'
#' @name trip-methods
#' @aliases trip-methods trip trip,SpatialPointsDataFrame,ANY-method
#' trip,ANY,TimeOrderedRecords-method trip,trip,ANY-method
#' trip,trip,TimeOrderedRecords-method [,trip-method [,trip,ANY,ANY,ANY-method [[<-,trip,ANY,missing-method
#' @param obj A \code{\link[sp]{SpatialPointsDataFrame}}, or an object that can
#' be coerced to one, containing at least two columns with the DateTime and ID
#' data as per \code{TORnames}.  It can also be a \code{trip} object for
#' redefining \code{TORnames}.
#' @param TORnames Either a \code{TimeOrderedRecords} object, or a 2-element
#' character vector specifying the DateTime and ID column of \code{obj}
#' @return
#'
#' A trip object, with the usual slots of a
#' \code{\link[sp]{SpatialPointsDataFrame}} and the added
#' \code{TimeOrderedRecords}. For the most part this can be treated as a
#' \code{data.frame} with \code{Spatial} coordinates.
#' @section Methods:
#'
#' Most of the methods available are by virtue of the sp package.  Some, such
#' as \code{split.data.frame} have been added to SPDF so that trip has the same
#' functionality.
#'
#' \describe{
#'
#' \item{trip}{\code{signature(obj="SpatialPointsDataFrame",
#' TORnames="ANY")}}The main construction.
#'
#' \item{trip}{\code{signature(obj="ANY", TORnames="TimeOrderedRecords")}:
#' create a \code{trip} object from a data frame.}
#'
#' \item{trip}{\code{signature(obj="trip", TORnames="ANY")}: (Re)-create a
#' \code{trip} object using a character vector for \code{TORnames}.}
#'
#' \item{trip}{\code{signature(obj="trip", TORnames="TimeOrderedRecords")}:
#' (re)-create a trip object using a \code{TimeOrderedRecords} object.}
#'
#' }
#' @seealso
#'
#' \code{\link{speedfilter}}, and \code{\link{tripGrid}} for simple(istic)
#' speed filtering and spatial time spent gridding.
#' @export
#' @examples
#'
#'
#' d <- data.frame(x=1:10, y=rnorm(10), tms=Sys.time() + 1:10, id=gl(2, 5))
#' coordinates(d) <- ~x+y
#' ## this avoids complaints later, but these are not real track data (!)
#' proj4string(d) <- CRS("+proj=laea +ellps=sphere")
#' (tr <- trip(d, c("tms", "id")))
#'
#' ## don't want adehabitatMA to be loaded as a requirement here
#' \dontrun{
#' ## a simple example with the common fixes required for basic track data
#'
#' dat <- read.csv("trackfile.csv")
#' names(dat)  ## e.g. [1] "long" "lat" "seal" "date" "local" "lq"
#' library(sp)
#' coordinates(dat) <- c("long", "lat")
#'
#' ## date/times may be in a particular time zone, please check
#' dat$gmt <- as.POSIXct(strptime(paste(dat$date, dat$local),
#'                       "%d-%b-%y %H:%M:%S"), tz="GMT")
#'
#' ## if there are problems in the data, this will error
#' tr <- trip(dat, c("gmt", "seal"))
#'
#' ## the following code tries to fix common problems
#'
#' ## remove completely-duplicated rows
#' dat <- dat[!duplicated(dat), ]
#' ## order the rows by seal, then by time
#' dat <- dat[order(dat$seal, dat$gmt), ]
#' ## fudge duplicated times
#' dat$gmt <- adjust.duplicateTimes(dat$gmt, dat$seal)
#'
#' ## finally, convert to Spatial and create trip object
#' coordinates(dat) <- c("long", "lat")
#' tr <- trip(dat, c("gmt", "seal"))
#' }
#'
#'
#' \dontrun{
#'    if (require(adehabitatLT)) {
#'      data(porpoise)
#'      porpoise <- as.trip(porpoise)
#'      proj4string(porpoise) <- CRS("+proj=utm +zone=21 +ellps=WGS84 +units=m +no_defs")
#'      summary(porpoise)
#'
#'    }
#'
#'
#'    ## extended example to check that our projection metadata is correct
#'    library(maptools)
#'    data(wrld_simpl)
#'    library(rgeos)
#'    library(raster)
#'
#'    ## 3 degrees either side (for half a zone . . .)
#'    ext <- as(extent(spTransform(porpoise, CRS(proj4string(wrld_simpl)))) + 3, "SpatialPolygons")
#'    proj4string(ext) <- CRS(proj4string(wrld_simpl))
#'    ## crop to the buffered tracks, and project to its native CRS
#'    w <- spTransform(gIntersection(wrld_simpl[grep("United States", wrld_simpl$NAME), ], ext),
#'     CRS(proj4string(porpoise)))
#'
#'    plot(w)
#'    lines(porpoise)
#' }
setGeneric("trip",
             function(obj, TORnames) standardGeneric("trip"))

if (!isGeneric("points"))
  setGeneric("points",
             function(x, ...) standardGeneric("points"))

if (!isGeneric("lines"))
  setGeneric("lines",
             function(x, ...) standardGeneric("lines"))

if (!isGeneric("text"))
  setGeneric("text",
             function(x, ...) standardGeneric("text"))

if (!isGeneric("subset"))
  setGeneric("subset",
             function(x, ...) standardGeneric("subset"))




##' TimeOrderedRecords
##'
##' Object to identify DateTimes and IDs in a Spatial object.
##'
##' @param x Character vector of 2 elements specifying the data columns of DateTimes and IDs
##' @return  \code{TimeOrderedRecords} holds a 2-element character vector, naming the data columns
##' of DateTimes and IDs.
##' @export
##' @examples
##' ##' tor <- TimeOrderedRecords(c("datetime", "ID"))
TimeOrderedRecords <- function(x) {
    new("TimeOrderedRecords", TOR.columns=x)
}



#'
#' Functions to retrieve DateTime and ID data from within (Spatial) data
#' frames.
#'
#'
#' Functions for retrieving the names of the columns used for DateTime and ID,
#' as well as the data.
#'
#' @name trip-accessors
#' @aliases trip-accessors getTORnames getTimeID
#' @param obj \code{trip} object.
#' @return
#'
#' \code{getTORnames} retrieves the column names from an object extending the
#' class \code{TimeOrderedRecords}, and \code{getTimeID} returns the data as a
#' data frame from an object extending the class \code{TimeOrderedRecords}.
#' @seealso
#'
#' \code{\link{trip-class}}, for the use of this class with
#' \code{\link[sp]{SpatialPointsDataFrame}}.
#'
#' \code{\link{trip}}
#' @keywords manip
#' @examples
#'
#'
#' tor <- TimeOrderedRecords(c("time", "id"))
#' getTORnames(tor)
#'
NULL

#' @rdname trip-accessors
#' @export
getTORnames <- function(obj) obj@TOR.columns

##' @rdname trip-accessors
##' @export
getTimeID <- function(obj) as.data.frame(obj)[, getTORnames(obj)]


setMethod("trip", signature(obj="SpatialPointsDataFrame", TORnames="ANY"),
          function(obj, TORnames) {
              if (is.factor(obj[[TORnames[2]]]))
                  obj[[TORnames[2]]] <- factor(obj[[TORnames[2]]])
              new("trip", obj, TimeOrderedRecords(TORnames))
          })

setMethod("trip", signature(obj="ANY", TORnames="TimeOrderedRecords"),
          function(obj, TORnames) {
              new("trip", obj, TORnames)
          })

setMethod("trip", signature(obj="trip", TORnames="TimeOrderedRecords"),
          function(obj, TORnames) {
              new("trip",
                  as(obj, "SpatialPointsDataFrame"),
                  TORnames)
          })

setMethod("trip", signature(obj="trip", TORnames="ANY"),
          function(obj, TORnames) {
              trip(as(obj, "SpatialPointsDataFrame"), TORnames)
          })

setReplaceMethod("[[",
                 signature(x="trip", i="ANY", j="missing", value="ANY"),
                 function(x, i, j, value) {
                     tor <- getTORnames(x)
                     x <- as(x, "SpatialPointsDataFrame")
                     x[[i]] <- value
                     trip(x, tor)
                 })

## S3 versions
dim.trip <- function(x) dim(as(x, "SpatialPointsDataFrame"))

as.data.frame.trip <- function(x, ...) {
    as.data.frame(as(x, "SpatialPointsDataFrame"), ...)
}

names.trip <- function(x) names(as(x, "SpatialPointsDataFrame"))

"names<-.trip" <- function(x, value) {
    names(x@data) <- value
    x@TOR.columns <- value
    x
}


###_ + sp methods

setMethod("points", signature(x="trip"),
          function(x, ...) points(as(x, "SpatialPointsDataFrame"), ...))
setMethod("text", signature(x="trip"),
          function(x, ...) text(as(x, "SpatialPointsDataFrame"), ...))

#setMethod("split", "SpatialPointsDataFrame", split.data.frame)

## setMethod("spTransform", signature=signature(x="trip", CRSobj="CRS"),
##           function(x, CRSobj, ...) tripTransform(x, CRSobj, ...))

## setMethod("spTransform", signature=signature(x="trip", CRSobj="character"),
##           function(x, CRSobj, ...) tripTransform(x, CRSobj, ...))

#' @exportMethod lines
setMethod("lines", signature(x="trip"),
          function(x,
                   col=hsv(seq(0, 0.9, length = length(unique(x[[getTORnames(x)[2]]]))),
                     0.8, 0.95),
                   ...) {
              plot(as(x, "SpatialLinesDataFrame"),  col=col, add=TRUE, ...)

          })
#' @exportMethod  plot
setMethod("plot", signature(x="trip", y="missing"),
          function(x, y, ...) {
              plot(as(x, "SpatialPoints"), ...)
          })


###_ + Subsetting trip

#' @exportMethod subset
setMethod("subset", signature(x="trip"),
          function(x,  ...) {
              spdf <- subset(as(x, "SpatialPointsDataFrame"), ...)
              tor <- getTORnames(x)
              if ( is.factor(spdf[[tor[2]]]))
                  spdf[[tor[2]]] <- factor(spdf[[tor[2]]])
              if (any(is.na(match(tor, names(spdf))))) {
                  msg <- paste("trip-defining Date or ID columns dropped,",
                               "reverting to SpatialPointsDataFrame\n\n")
                  cat(msg)
                  return(spdf)
              } else {
                  tst <- any(tapply(spdf[[tor[1]]],
                                    spdf[[tor[2]]], length) < 3)
                  if (tst) {
                      msg <- paste("subset loses too many locations,",
                               "reverting to SpatialPointsDataFrame\n\n")
                      cat(msg)
                      return(spdf)
                  } else return(trip(spdf, tor))
              }
          })

##' @param x trip object
##' @param i,j,\dots indices specifying elements to extract 
##' @param drop unused but necessary for method consistency
##' @rdname trip-methods
setMethod("[", signature(x="trip"),
          function(x, i, j, ... , drop=TRUE) {
              missing.i <- missing(i)
              missing.j <- missing(j)
              nargs <- nargs() # e.g., a[3,] gives 2 for nargs, a[3] gives 1.
              if (missing.i && missing.j) {
                  i <- j <- TRUE
              } else if (missing.j && !missing.i) {
                  if (nargs == 2) {
                      j <- i; i <- TRUE
                  } else j <- TRUE
              } else if (missing.i && !missing.j) i <- TRUE
              if (is.matrix(i)) {
                  msg <- paste("matrix argument not supported in",
                               "SpatialPointsDataFrame selection")
                  stop(msg)
              }
              if (any(is.na(i)))
                  stop("NAs not permitted in row index")
              spdf <- as(x, "SpatialPointsDataFrame")[i, j, ..., drop=drop]
              tor <- getTORnames(x)
              if (is.factor(spdf[[tor[2]]]))
                  spdf[[tor[2]]] <- factor(spdf[[tor[2]]])
              if (any(is.na(match(tor, names(spdf))))) {
                  msg <- paste("trip-defining Date or ID columns dropped,",
                               "reverting to SpatialPointsDataFrame\n\n")
                  cat(msg)
                  return(spdf)
              } else {
                  tst <- any(tapply(spdf[[tor[1]]],
                                    spdf[[tor[2]]], length) < 3)
                  if (tst) {
                      msg <- paste("subset loses too many locations,",
                                   "reverting to SpatialPointsDataFrame\n\n")
                      cat(msg)
                      return(spdf)
                  } else {
                      return(trip(spdf, tor))
                  }
              }
          })


###_ + Summary, print, and show

#' @exportMethod summary
setMethod("summary", signature(object="trip"),
          function(object, ...) {
              obj <- list(spdf=summary(as(object,
                            "SpatialPointsDataFrame")))
              tids <- getTimeID(object)
              time <- tids[, 1]
              ids <- tids[, 2]
              ## list of distances only, km/hr or units of projection
              dists <- .distances(object)
              rmsspeed <- split(speedfilter(object, max.speed = 1, test = TRUE)$rms, ids)

              ## list of time diferences only, in hours
              dtimes <- lapply(split(time, ids), function(x) diff(unclass(x)/3600))
              speeds <- vector("list", length(dtimes))
              for (i in seq_along(speeds)) speeds[[i]] <- dists[[i]] / dtimes[[i]]

              obj <- within(obj, {
                  class <- class(object)
                  tmins <- tapply(time, ids, min) +
                      ISOdatetime(1970, 1, 1, 0, 0,0, tz="GMT")
                  tmaxs <- tapply(time, ids, max) +
                      ISOdatetime(1970, 1, 1, 0, 0,0, tz="GMT")
                  tripID <- levels(factor(ids))
                  nRecords <- tapply(time, ids, length)
                  TORnames <- getTORnames(object)
                  tripDuration <- tapply(time, ids, function(x) {
                      x <- format(diff(range(x)))
                  })
                  tripDurationSeconds <- tapply(time, ids, function(x) {
                      x <- diff(range(unclass(x)))
                  }
                                                )
                  tripDistance <- sapply(dists, sum)
                  meanSpeed <- sapply(speeds, mean)
                  maxSpeed <- sapply(speeds, max)
                  meanRMSspeed <- sapply(rmsspeed, mean, na.rm = TRUE)
                  maxRMSspeed <- sapply(rmsspeed, max, na.rm = TRUE)
              })
              class(obj) <- "summary.TORdata"
              ## invisible(obj)
              obj
          })

as.data.frame.summary.TORdata <- function(x, row.names = NULL, optional = FALSE, ...) {
        dsumm <- data.frame(tripID=x$tripID,
                        No.Records=x$nRecords,
                        startTime=x$tmins,
                        endTime=x$tmaxs,
                        tripDuration=x$tripDuration,
                        tripDistance=x$tripDistance,
                        meanSpeed = x$meanSpeed,
                        maxSpeed = x$maxSpeed,
                        meanRMSspeed = x$meanRMSspeed,
                        maxRMSspeed = x$maxRMSspeed)
  dsumm
}

#' @rdname trip-accessors
#' @method print summary.TORdata
#' @param x trip object
#' @param \dots currently ignored
#' @export
print.summary.TORdata <- function(x, ...) {
    dsumm <- as.data.frame(x)
  torns <- x[["TORnames"]]
    names(dsumm)[1] <- paste(names(dsumm)[1],
                             " (\"", torns[2], "\")", sep="")
    names(dsumm)[3] <- paste(names(dsumm)[3],
                             " (\"", torns[1], "\")", sep="")
    names(dsumm)[4] <- paste(names(dsumm)[4],
                             " (\"", torns[1], "\")", sep="")


    rownames(dsumm) <- seq(nrow(dsumm))
    ## dsumm <- as.data.frame(lapply(dsumm, as.character))
    cat(paste("\nObject of class ", x[["class"]], "\n", sep=""))
    print(format(dsumm, ...))
    tripDurationSeconds <- sum(x$tripDurationSeconds)
    tripDurationHours <- sum(x$tripDurationSeconds) / 3600
    cat(paste("\nTotal trip duration: ",
              tripDurationSeconds, " seconds (",
              as.integer(tripDurationHours), " hours, ",
              round((tripDurationHours -
                     as.integer(tripDurationHours)) * 3600),
              " seconds)\n", sep=""))
    cat(paste("\nDerived from Spatial data:\n\n", sep=""))
    print(x$spdf)
    cat("\n")
}

#' @exportMethod show
setMethod("show", signature(object="summary.TORdata"),
          function(object) print.summary.TORdata(object))

print.trip <- function(x, ...) {
    xs <- summary(x)
    dsumm <- data.frame(tripID=xs$tripID,
                        No.Records=xs$nRecords,
                        startTime=xs$tmins,
                        endTime=xs$tmaxs,
                        tripDuration=xs$tripDuration)
    torns <- xs[["TORnames"]]
    names(dsumm)[1] <- paste(names(dsumm)[1], " (\"",
                             torns[2], "\")", sep="")
    names(dsumm)[3] <- paste(names(dsumm)[3], " (\"",
                             torns[1], "\")", sep="")
    names(dsumm)[4] <- paste(names(dsumm)[4], " (\"",
                             torns[1], "\")", sep="")
    rownames(dsumm) <- 1:nrow(dsumm)
    ## dsumm <- as.data.frame(lapply(dsumm, as.character))
    cat(paste("\nObject of class ", xs[["class"]], "\n", sep=""))
    print(format(dsumm, ...))
    cat("\n")
    nms <- names(x)
    clss <- unlist(lapply(as.data.frame(x@data), function(x) class(x)[1]))
    sdf <- data.frame(data.columns=nms, data.class=clss)
    sdf[[" "]] <- rep("", nrow(sdf))
    sdf[[" "]][nms == torns[1]] <- "**trip DateTime**"
    sdf[[" "]][nms == torns[2]] <- "**trip ID**      "
    row.names(sdf) <- seq(nrow(sdf))
    print(sdf)
    cat("\n")
}

setMethod("show", signature(object="trip"),
          function(object) print.trip(object))

setMethod("recenter", signature(obj="trip"),
          function(obj) {
              proj <- is.projected(obj)
              if (is.na(proj)) {
                  msg <- paste("unknown coordinate reference system:",
                               "assuming longlat")
                  warning(msg)
                  ## projargs <- CRS("+proj=longlat")
              }
              if (!is.na(proj) & proj) {
                  msg <- paste("cannot recenter projected coordinate",
                               "reference system")
                  stop(msg)
              }
              projargs <- CRS(proj4string(obj))
              crds <- coordinates(obj)
              inout <- (crds[, 1] < 0)
              if (all(inout)) {
                  crds[, 1] <- crds[, 1] + 360
                  if (!is.na(proj)) projargs <- CRS(paste(proj4string(obj),
                                                          "+over"))
              } else {
                  if (any(inout)) {
                      crds[, 1] <- ifelse(inout, crds[, 1] + 360,
                                          crds[, 1])
                      if (!is.na(proj))
                          projargs <- CRS(paste(proj4string(obj), "+over"))
                  }
              }
              trip(new("SpatialPointsDataFrame",
                       SpatialPoints(crds, projargs),
                       data=obj@data, coords.nrs=obj@coords.nrs),
                   obj@TOR.columns)
          })




setMethod("spTransform", signature("trip", "CRS"),
          function(x, CRSobj, ...) {
            if (!("rgdal" %in% loadedNamespaces())) {
              ns <- try(loadNamespace("rgdal"))
              if (isNamespace(ns)) {
                message("[loaded the rgdal namespace]")
              } else {
                msg <- paste("This method requires the rgdal package",
                             "but is unable to load rgdal namespace",
                             sep=",")
                stop(msg)
              }
            }
            pts <- spTransform(as(x, "SpatialPointsDataFrame"),
                               CRSobj, ...)
            trip(pts, getTORnames(x))
          })

## method to allow transformation with character only
setMethod("spTransform", signature("Spatial", "character"), 
          function(x, CRSobj, ...) {
            
            .local <- function (object, pstring, ...) 
            {
              crs <- try(CRS(pstring))
              if (inherits(crs, "try-error")) { stop(sprintf("cannot determine valid CRS from %s", pstring))
              } else {
                spTransform(x, crs)
              }
            }
            
            .local(x, pstring = CRSobj, ...)
            
          })

