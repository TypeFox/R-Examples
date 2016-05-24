
#' Function to ensure dates and times are in order with trip ID
#'
#'
#' A convenience function, that removes duplicate rows, sorts by the date-times
#' within ID, and removes duplicates from a data frame or
#' SpatialPointsDataFrame.
#'
#'
#' @param x \code{\link{data.frame}} or
#' \code{\link[sp]{SpatialPointsDataFrame}}
#' @param tor character vector of names of date-times and trip ID columns
#' @return \code{\link{data.frame}} or
#' \code{\link[sp]{SpatialPointsDataFrame}}.
#' @note
#'
#' It's really important that data used are of a given quality, but this
#' function makes the most common trip problems easy to apply.
#' @seealso \code{\link{trip}}
#' @export forceCompliance
forceCompliance <- function(x, tor) {
    isSpatial <- is(x, "SpatialPointsDataFrame")
    if (isSpatial) {
        crd.nrs <- x@coords.nrs
        x <- as.data.frame(x)
    }
    levs <- unique(x[[tor[2]]])
    tooshort <- tapply(x[[1]], x[[tor[2]]], function(x) length(x) < 3)
    x <- x[x[[tor[2]]] %in% levs[!tooshort], ]
    x <- x[!duplicated(x), ]
    x <- x[order(x[[tor[2]]], x[[tor[1]]]), ]
    x[[tor[1]]] <- adjust.duplicateTimes(x[[tor[1]]], x[[tor[2]]])
    if (isSpatial) {
        coordinates(x) <- crd.nrs
        x <- trip(x, tor)
    }
    x
}

##' @rdname trip-internal
.intpFun <- function(x) {
    len <- round(x[3] + 1)
    new <- seq(x[1], x[2], length=len)
    if (len > 1)
        new[-len]
    else new
}

interpequal <- function(x, dur=NULL, quiet=FALSE) {
  if (!is(x, "trip"))
    stop("only trip objects supported")
  if (is.null(dur))
    stop("equal time duration must be specified \"dur=?\"")
  ## x must be a single trip
  tor <- getTORnames(x)
  tids <- getTimeID(x)
  time <- tids[, 1]
  id <- factor(tids[, 2])
  coords <- coordinates(x)
  x <- coords[,1]
  y <- coords[,2]
  levs <- levels(id)
  newPts <- vector("list", length(levs))
  ##if (is.null(dur))
  ##   dur <- as.numeric(min(unlist(tapply(as.integer(time),
  ##            id, diff))))
  for (isub in seq_along(levs)) {
    ind <- id == levs[isub]
    xx <- x[ind]
    yy <- y[ind]
    tms <- time[ind]
    dt <- diff(as.numeric(tms))
    dtn <- dt/dur
    ax <- cbind(xx, c(xx[-1], xx[length(xx)]), c(dtn, 0))
    ay <- cbind(yy, c(yy[-1], yy[length(yy)]), c(dtn, 0))
    intime <- as.numeric(tms) - min(as.numeric(tms))
    at <- cbind(intime, c(intime[-1], intime[length(intime)]),
                c(dtn, 0))
    nx <- unlist(apply(ax, 1, .intpFun))
    ny <- unlist(apply(ay, 1, .intpFun))
    nt <- unlist(apply(at, 1, .intpFun)) + min(tms)
    ni <- factor(rep(levs[isub], length=length(nt)))
##    newPts <- rbind(newPts,
    newPts[[isub]] <- data.frame(x=nx, y=ny, time=nt, id=ni)
  }
  newPts <- do.call("rbind", newPts)

  origTotal <- sum(tapply(time, id, function(x) {
    diff(range(as.numeric(x)))
  }))
  newTotal <- nrow(newPts) * dur
  uType <- "hours"
  hTotal <- sum(tapply(time, id, function(x) {
    difftime(range(x)[2], range(x)[1], units=uType)
  }))
  if (!quiet) {
    cat("lost seconds=", as.integer(origTotal - newTotal),
        " out of a total ", hTotal, " ", uType, "\n")
  }
  coordinates(newPts) <- c("x", "y")
  names(newPts) <- tor
  newPts
}






#' Adjust duplicate DateTime values
#'
#'
#' Duplicated DateTime values within ID are adjusted forward (recursively) by
#' one second until no duplicates are present. This is considered reasonable
#' way of avoiding the nonsensical problem of duplicate times.
#'
#'
#' This function is used to remove duplicate time records in animal track data,
#' rather than removing the record completely.
#'
#' @param time vector of DateTime values
#' @param id vector of ID values, matching DateTimes that are assumed sorted
#' within ID
#' @return
#'
#' The adjusted DateTime vector is returned.
#' @section Warning:
#'
#' I have no idea what goes on at CLS when they output data that are either not
#' ordered by time or have duplicates. If this problem exists in your data it's
#' probably worth finding out why.
#' @seealso \code{\link{readArgos}}
#' @examples
#'
#'
#' ## DateTimes with a duplicate within ID
#' tms <- Sys.time() + c(1:6, 6, 7:10) *10
#' id <- rep("a", length(tms))
#' range(diff(tms))
#'
#' ## duplicate record is now moved one second forward
#' tms.adj <- adjust.duplicateTimes(tms, id)
#' range(diff(tms.adj))
#'
#'
#' @export adjust.duplicateTimes
adjust.duplicateTimes <- function (time, id) {
    dups <- unlist(tapply(time, id, duplicated), use.names=FALSE)
    if (any(dups)) {
        time[dups] <- time[dups] + 1
        time <- Recall(time, id)
    }
    time
}



#' Assign numeric values for Argos "class"
#'
#'
#' Assign numeric values for Argos "class" by matching the levels available to
#' given numbers. An adjustment is made to allow sigma to be specified in
#' kilometeres, and the values returned are the approximate values for longlat
#' degrees.  It is assumed that the levels are part of an "ordered" factor from
#' least precise to most precise.
#'
#'
#' The available levels in Argos are \code{levels=c("Z", "B", "A", "0", "1",
#' "2", "3")}.
#'
#' The actual sigma values given by default are (as far as can be determined) a
#' reasonable stab at what Argos believes.
#'
#' @param x factor of Argos location quality "classes"
#' @param sigma numeric values (by default in kilometres)
#' @param adjust a numeric adjustment to convert from kms to degrees
#' @return
#'
#' Numeric values for given levels.
#' @keywords manip
#' @examples
#'
#'
#' cls <- ordered(sample(c("Z", "B", "A", "0", "1", "2", "3"), 30,
#'                       replace=TRUE),
#'                levels=c("Z", "B", "A", "0", "1", "2", "3"))
#' argos.sigma(cls)
#'
#'
#' @export argos.sigma
argos.sigma <- function(x, sigma=c(100, 80, 50, 20, 10, 4,  2),
                        adjust=111.12) {
    sigma <- sigma / adjust
    names(sigma) <- levels(x)
    sigma[x]
}



#' Read Argos "DAT" or "DIAG" files
#'
#'
#' Return a (Spatial) data frame of location records from raw Argos files.
#' Multiple files may be read, and each set of records is appended to the data
#' frame in turn.  Basic validation of the data is enforced by default.
#'
#'
#' \code{readArgos} performs basic validation checks for class \code{trip} are
#' made, and enforced based on \code{correct.all}:
#'
#' No duplicate records in the data, these are simply removed.  Records are
#' ordered by DateTime ("date", "time", "gmt") within ID ("ptt").  No duplicate
#' DateTime values within ID are allowed: to enforce this the time values are
#' moved forward by one second - this is done recursively and is not robust.
#'
#' If validation fails the function will return a
#' \code{\link[sp]{SpatialPointsDataFrame}}.  Files that are not obviously of
#' the required format are skipped.
#'
#' Argos location quality data "class" are ordered, assuming that the available
#' levels is \code{levels=c("Z", "B", "A", "0", "1", "2", "3")}.
#'
#' A projection string is added to the data, assuming the PROJ.4 longlat - if
#' any longitudes are greater than 360 the PROJ.4 argument "+over" is added.
#'
#' \code{readDiag} simply builds a \code{data.frame}.
#'
#' @aliases readArgos readDiag
#' @param x vector of file names of Argos "DAT" or "DIAG" files.
#' @param correct.all logical - enforce validity of data as much as possible?
#' (see Details)
#' @param dtFormat the DateTime format used by the Argos data "date" and "time"
#' pasted together
#' @param tz timezone - GMT/UTC is assumed
#' @param duplicateTimes.eps what is the tolerance for times being duplicate?
#' @param p4 PROJ.4 projection string, "+proj=longlat +ellps=WGS84" is assumed
#' @param verbose if TRUE, details on date-time adjustment is reported
#' @return
#'
#' \code{readArgos} returns a \code{trip} object, if all goes well, or simply a
#' \code{\link[sp]{SpatialPointsDataFrame}}.
#'
#' \code{readDiag} returns a \code{data.frame} with 8 columns:
#' \itemize{
#' \item {\code{lon1},\code{lat1} first pair of coordinates}
#' \item {\code{lon1},\code{lat1} second pair of coordinates}
#' \item {gmt DateTimes as POSIXct}
#' \item {id Platform Transmitting Terminal (PTT) ID}
#' \item {lq Argos location quality class}
#' \item {iq some other thing}
#' }
#' @section Warning :
#'
#' This works on some Argos files I have seen, it is not a guaranteed method
#' and is in no way linked officially to Argos.
#' @seealso
#'
#' \code{\link{trip}}, \code{\link[sp]{SpatialPointsDataFrame}},
#' \code{\link{adjust.duplicateTimes}}, for manipulating these data, and
#' \code{\link{argos.sigma}} for relating a numeric value to Argos quality
#' "classes".
#'
#' \code{\link{sepIdGaps}} for splitting the IDs in these data on some minimum
#' gap.
#'
#' \code{\link{order}}, \code{\link{duplicated}}, , \code{\link{ordered}} for
#' general manipulation of this type.
#' @references
#'
#' The Argos data documentation was (ca. 2003) at
#' http://www.argos-system.org/manual.  Specific details on the PRV
#' ("provide data") format were found in Chapter 4_4_8, originally at 
#' 'http://www.cls.fr/manuel/html/chap4/chap4_4_8.htm'.
#' @keywords IO manip
#' @export readArgos
readArgos <- function (x, correct.all=TRUE, dtFormat="%Y-%m-%d %H:%M:%S",
                       tz="GMT", duplicateTimes.eps=1e-2,
                       p4="+proj=longlat +ellps=WGS84", verbose=FALSE) {
    ## add "correct.all" argument - just return data frame if it fails, with
    ## suggestions of how to sort/fix it


  ## this should be heaps faster
    dout <- vector("list", length(x))
    for (icon in seq_along(x)) {
        old.opt <- options(warn=-1)
        dlines <- strsplit(readLines(x[icon]), "\\s+", perl=TRUE)
        options(old.opt)
        loclines <- sapply(dlines, length) == 12
        if (any(loclines)) {
            dfm <- matrix(unlist(dlines[sapply(dlines, length) == 12]),
                          ncol=12, byrow=TRUE)
            if (dfm[1,7] == "LC") {
                msg <- paste(" appears to be a diag file, skipping.",
                             "Use readDiag to obtain a dataframe. \n\n")
            	cat("file ", icon, msg)
            	next
            }
            df <- vector("list", 12)
            names(df) <- c("prognum", "ptt", "nlines", "nsensor",
                           "satname", "class", "date", "time", "latitude",
                           "longitude", "altitude", "transfreq")
            for (i in c(1:4, 9:12)) df[[i]] <- as.numeric(dfm[, i])
            for (i in 5:6) df[[i]] <- factor(dfm[, i])
            for (i in 7:8) df[[i]] <- dfm[, i]
            df <- as.data.frame(df)
            df$gmt <- as.POSIXct(strptime(paste(df$date, df$time),
                                          dtFormat), tz)
            dout[[icon]] <- df
        } else {
            cat("Problem with file: ", x[icon], " skipping\n")

        }
    }
    if (all(sapply(dout, is.null)))
        stop("No data to return: check the files")

    dout <- do.call(rbind, dout)
    if (correct.all) {
        ## should add a reporting mechanism for these as well
        ##  and return a data.frame if any of the tests fail
        ## sort them
        dout <- dout[order(dout$ptt, dout$gmt), ]
        ## remove duplicate rows
        dout <- dout[!duplicated(dout), ]
        ## adjust duplicate times (now that they are sorted properly)
        dt.by.id <- unlist(tapply(dout$gmt, dout$ptt,
                                  function(x) c(-1, diff(x))))
        dup.by.eps <- which(abs(dt.by.id) < duplicateTimes.eps)
        if (length(dup.by.eps) >= 1) {
            if (verbose) {
                cat("Adjusting duplicate times\n.....\n")
                for (i in  dup.by.eps) {
                    ind <- i + (-2:1)
                    print(cbind(dout[ind,c("ptt", "gmt", "class")],
                                row.number=ind))
                }
            }
            dout$gmt <- adjust.duplicateTimes(dout$gmt, dout$ptt)
            if (verbose) {
                cat("\n  Adjusted records now: \n\n")
                for (i in  dup.by.eps) {
                    ind <- i + (-2:1)
                    print(cbind(dout[ind,c("ptt", "gmt", "class")],
                                row.number=ind))
                }
            }
        }
        if(any(dout$longitude > 180)) {
            msg <- paste("\nLongitudes contain values greater than 180,",
                         "assuming proj.4 +over\n\n")
            cat(msg)
            p4 <- "+proj=longlat +ellps=WGS84 +over"
        }
        dout$class <- ordered(dout$class,
                              levels=c("Z", "B", "A", "0", "1", "2", "3"))
        coordinates(dout) <- c("longitude", "latitude")
        proj4string(dout) <- CRS(p4)
        ##tor <- TimeOrderedRecords(c("gmt", "ptt"))
        test <- try(dout <- trip(dout, c("gmt", "ptt")))
        if (!is(test, "trip")) {
            cat("\n\n\n Data not validated: returning object of class ",
                class(dout), "\n")
            return(dout)
        }
        ## for now, only return spdftor if correct.all is TRUE
        cat("\n\n\n Data fully validated: returning object of class ",
            class(dout), "\n")
        return(dout)
    }
    cat("\n\n\n Data not validated: returning object of class ",
        class(dout), "\n")
    dout
}

##' @rdname readArgos
##' @export
readDiag <- function (x) {
  data <- NULL
  for (fl in x) {
    d <- readLines(fl)
    locs <- d[grep("LON1", d, ignore.case=TRUE)]
    tms <- d[grep("DATE", d, ignore.case=TRUE)]
    bad <- (grep("\\?", locs))
    if (length(bad) > 0) {
      if (length(locs[-bad]) == 0) {
        warning(paste("no valid locations in:", fl, "\n ...ignoring"))
        next
      }
      locs <- locs[-bad]
      tms <- tms[-(bad)]
    }
    dlines <- paste(locs, tms)
    dlines <- strsplit(dlines, "\\s+", perl=TRUE)
    reclen <- length(dlines[[1]])
    dfm <- matrix(unlist(dlines[sapply(dlines, length) ==
                                  reclen]), ncol=reclen, byrow=TRUE)
    lonlat <- dfm[, c(4, 7, 10, 13)]
    dic <- dfm[, c(14, 17, 18, 21, 24), drop=FALSE]
    id <- dic[, 1]
    gmt <- as.POSIXct(strptime(paste(dic[, 2], dic[, 3]),
                               "%d.%m.%y %H:%M:%S"), tz="GMT")
    lq <- dic[, 4]
    iq <- dic[, 5]
    ll <- as.vector(lonlat)
    ll[grep("S", ll)] <- paste("-", ll[grep("S", ll)], sep="")
    ll <- gsub("S", "", ll)
    ll[grep("N", ll)] <- paste("", ll[grep("N", ll)], sep="")
    ll <- gsub("N", "", ll)
    ll[grep("E", ll)] <- paste("", ll[grep("E", ll)], sep="")
    ll <- gsub("E", "", ll)
    ll[grep("W", ll)] <- paste("-", ll[grep("W", ll)], sep="")
    ll <- gsub("W", "", ll)
    ll <- matrix(as.numeric(ll), ncol=4)
    lon <- ll[, 2]
    lon2 <- ll[, 4]
    lq <- factor(lq, ordered=TRUE,
                 levels=c("Z", "B", "A", "0", "1", "2", "3"))
    data <- rbind(data,
                  data.frame(lon1=lon, lat1=ll[, 1],
                             lon2=lon2, lat2=ll[, 3],
                             gmt=gmt, id=id, lq=lq, iq=iq))
  }
  data
}



#' Separate a set of IDs based on gaps
#'
#'
#' A new set of ID levels can be created by separating those given based on a
#' minimum gap in another set of data. This is useful for separating
#' instruments identified only by their ID into separate events in time.
#'
#'
#' The assumption is that a week is a long time for a tag not to record
#' anything.
#'
#' @param id existing ID levels
#' @param gapdata data matching \code{id} with gaps to use as separators
#' @param minGap the minimum "gap" to use in gapdata to create a new ID level
#' @return
#'
#' A new set of ID levels, named following the pattern that "ID" split into 3
#' would provided "ID", "ID\_2" and "ID\_3".
#' @section Warning:
#'
#' It is assumed that each vector provides is sorted by \code{gapdata} within
#' \code{id}. No checking is done, and so it is suggested that this only be
#' used on ID columns within existing, validated \code{trip} objects.
#' @seealso \code{\link{trip}}
#' @keywords manip
#' @examples
#'
#'
#' id <- gl(2, 8)
#' gd <- Sys.time() + 1:16
#' gd[c(4:6, 12:16)] <- gd[c(4:6, 12:16)] + 10000
#' sepIdGaps(id, gd, 1000)
#'
#'
#' @export sepIdGaps
sepIdGaps <- function(id, gapdata, minGap=3600 * 24 * 7) {
    toSep <- tapply(gapdata, id,
                    function(x) which(diff(unclass(x) ) > minGap))
    tripID <- split(as.character(id), id)
    for (i in 1:length(tripID)) {
        this <- toSep[[i]]
        thisID <- tripID[[i]][1]
        if (length(this) > 0) {
            for (n in 1:length(this)) {
                tripID[[i]][(this[n]+1):length(tripID[[i]])] <-
                    paste(thisID, n + 1, sep="_")
            }
        }
    }
    unsplit(tripID, id)
}


## taken from package sp/src/gcdist.c

##' @rdname trip-internal
.distances <- function(x) {
  proj <- is.projected(x)
  if (is.na(proj)) proj <- FALSE


  lapply(split(x, x[[getTORnames(x)[2]]]), function(x) trackDistance(coordinates(x), longlat = !proj))

}

##' @rdname trip-internal
.gcdist.c <- function(lon1, lat1, lon2, lat2) {
    DE2RA <- pi / 180
    a <- 6378.137            # /* WGS-84 equatorial radius in km */
    f <- 1.0 / 298.257223563 # /* WGS-84 ellipsoid flattening factor */
    lat1R <- lat1 * DE2RA
    lat2R <- lat2 * DE2RA
    lon1R <- lon1 * DE2RA
    lon2R <- lon2 * DE2RA
    F <- ( lat1R + lat2R ) / 2.0
    G <- ( lat1R - lat2R ) / 2.0
    L <- ( lon1R - lon2R ) / 2.0
    sinG2 <- sin( G ) ^ 2
    cosG2 <- cos( G ) ^ 2
    sinF2 <- sin( F ) ^ 2
    cosF2 <- cos( F ) ^ 2
    sinL2 <- sin( L ) ^ 2
    cosL2 <- cos( L ) ^ 2
    S <- sinG2 * cosL2 + cosF2 * sinL2
    C <- cosG2 * cosL2 + sinF2 * sinL2
    w <- atan( sqrt( S / C ) )
    R <- sqrt( S * C ) / w
    D <- 2 * w * a
    H1 <- ( 3 * R - 1 ) / ( 2 * C )
    H2 <- ( 3 * R + 2 ) / ( 2 * S )
    dist <- D * ( 1 + f * H1 * sinF2 * cosG2 - f * H2 * cosF2 * sinG2 )
    ##dist <- ifelse((abs(lat1 - lat2) < .Machine$double.eps) & (abs(lon1 - lon2) < .Machine$double.eps), 0.0, dist)
    #dist <- ifelse(abs((abs(lon1) + abs(lon2)) - 360.0) < .Machine$double.eps, 0.0, dist)
    dist[is.na(dist)] <- 0
    dist
}

##trackDistance <- function(x) UseMethod("trackDistance")


#' Determine distances along a track
#'
#'
#' Calculate the distances between subsequent 2-D coordinates using Euclidean
#' or Great Circle distance (WGS84 ellipsoid) methods.
#'
#'
#' If \code{x1} is a trip object, arguments \code{x2}, \code{x3}, \code{y2} are
#' ignored and the return result has an extra element for the start point of
#' each individual trip, with value 0.0.
#'
#' The \code{prev} argument is ignore unless x1 is a trip.
#'
#' Distance values are in the units of the input coordinate system when longlat
#' is FALSE, and in kilometres when longlat is TRUE.
#'
#' This originally used \code{\link[sp]{spDistsN1}} but now implements the sp
#' \code{gcdist} source directly in R.
#'
#' @aliases trackDistance trackDistance.default trackDistance.trip
#' @param x1 trip object, matrix of 2-columns, with x/y coordinates OR a vector
#' of x start coordinates
#' @param x2 vector of x end coordinates, if x1 is not a matrix
#' @param y1 vector of y start coordinates, if x1 is not a matrix
#' @param y2 vector of y end coordinates, if x1 is not a matrix
#' @param longlat if FALSE, Euclidean distance, if TRUE Great Circle distance
#' @param prev if TRUE and x1 is a trip, the return value has a padded end
#' value (\"prev\"ious), rather than start (\"next\")
#' @return Vector of distances between coordinates.
#' @references Original source taken from sp package.
#' @author Roger Bivand and Michael Sumner
#' @examples
#'   ## Continuing the example from '?"trip-methods"':
#' utils::example("trip-methods", package="trip",
#'                ask=FALSE, echo=FALSE)
#'
#'  ## the method knows this is a trip, so there is a distance for every
#'  ## point, including 0s as the start and at transitions between
#'  ## individual trips
#' trackDistance(tr)
#'
#' ## the default method does not know about the trips, so this is
#' ##(n-1) distances between all points
#' ## trackDistance(coordinates(tr), longlat = FALSE)
#'
#' ## we get NA at the start, end and at transitions between trips
#'
#'  \dontrun{
#'  require(rgdal)
#'  trackAngle(tr)
#'  }
#' @export trackDistance
trackDistance <- function(x1, y1, x2, y2, longlat=TRUE, prev = FALSE) UseMethod("trackDistance")

##' @export
trackDistance.default <- function(x1, y1, x2, y2, longlat=TRUE, prev = FALSE) {
    if (missing(y1)) {
        if (!is.matrix(x1))
            stop("x1 is not a matrix and multiple arguments not specified")
        if (nrow(x1) < 2) stop("x1 has too few rows")
        if (ncol(x1) < 2) stop("x1 has too few columns")
        x2 <- x1[-1, 1]
        y1 <- x1[-nrow(x1), 2]
        y2 <- x1[-1, 2]
        x1 <- x1[-nrow(x1), 1]
    }
    nx <- length(x1)
    if (nx != length(y1) | nx != length(x2) | nx != length(y2))
        stop("arguments must have equal lengths")
    if (longlat) {
        .gcdist.c(x1, y1, x2, y2)
    } else sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2)
}

##' @export
trackDistance.trip <- function(x1, y1, x2, y2, longlat = TRUE, prev = FALSE) {
    unlist(lapply(.distances(x1), function(x) if (prev) {c(x, 0)} else {c(0, x)}))
}




#' Determine internal angles along a track
#'
#'
#' Calculate the angles between subsequent 2-D coordinates using Great Circle
#' distance (spherical) methods.
#'
#' If \code{x} is a trip object, the return result has an extra element for the
#' start and end point of each individual trip, with value NA.
#'
#' This is an optimized hybrid of "raster::bearing" and
#' \code{\link[maptools]{gzAzimuth}}.
#'
#' @rdname trackAngle
#' @param x trip object, or matrix of 2-columns, with x/y coordinates
#' @return Vector of angles (degrees) between coordinates.
#' @rdname trackAngle
#' @export trackAngle
trackAngle <- function(x) {
  UseMethod("trackAngle")
}

#' @rdname trackAngle
#' @method trackAngle trip

#' @export
trackAngle.trip <- function(x) {
  isproj <- is.projected(x)
  if (is.na(isproj)) {
    warning("object CRS is NA, assuming longlat")
  } else {
    if (isproj) {
      x <- try(spTransform(x, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
      if (inherits(x, "try-error")) {
        stop("Object x is projected, attempts to transform to longlat failed. Is rgdal installed?")
      }
    }
  }
  st <- split(x, x[[getTORnames(x)[2]]])
  unlist(lapply(st, function(x1) c(NA, trackAngle(coordinates(x1)), NA)))

}

#' @rdname trackAngle
#' @method trackAngle default
#' @export
trackAngle.default <- function(x) {
  n <- nrow(x)
  ## MDSumner 2013-06-14 not sure what to expose here, will start with optimized gzAzimuth(abdali)/bearing() hybrid
  ##if (type == "geosphere") {
  ##  require(geosphere)
  ##  angle <- bearing(xy[2:(n-1),],xy[3:n,]) - bearing(xy[2:(n-1),],xy[1:(n-2),])
  ##} else {
  ##  if(!type == "abdali") stop(sprintf("type '%s' not implemented", type))
    angle <- .abdali(x[2:(n-1),],x[3:n,]) - .abdali(x[2:(n-1),],x[1:(n-2),])

  ##}
  angle <- (angle+180)%%360-180
  abs(angle)
}


## "abdali", replacement for raster::bearing
##' @rdname trip-internal
.abdali <- function (p1, p2)
{
  stopifnot(nrow(p1) == nrow(p2))

  toRad <- pi/180
  p1 <- p1 * toRad
  p2 <- p2 * toRad
  keep <- !(.rowSums(p1 == p2, nrow(p1), ncol(p1)) == 2L)
  res <- rep(as.numeric(NA), length = nrow(p1))
  if (sum(keep) == 0) {
    return(res)
  }
  p1 <- p1[keep, , drop = FALSE]
  p2 <- p2[keep, , drop = FALSE]
  dLon = p2[, 1] - p1[, 1]
  y = sin(dLon)
  x = cos(p1[, 2]) * tan(p2[, 2]) - sin(p1[, 2]) * cos(dLon)
  azm = atan2(y, x)/toRad
  res[keep] <- (azm + 360)%%360
  return(res)
}


##n <- 10000;x <- cbind(runif(n, -180, 180), runif(n, -90, 90));

##max(abs(trackAngle(x) - trackAngle(x, type = "abdali")))
## [1] 1.136868e-13

##library(rbenchmark)

##n <- 5000;x <- cbind(runif(n, -180, 180), runif(n, -90, 90));

##benchmark(trackAngle(x, type = "geosphere"), trackAngle(x, type = "abdali"), replications = 300)
##test replications elapsed relative user.self sys.self user.child sys.child
##2    trackAngle(x, type = "abdali")          300    1.62    1.000      1.62        0         NA        NA
##1 trackAngle(x, type = "geosphere")          300    8.49    5.241      8.49        0         NA        NA


## TODO:
## tidier!
##-----------------------------------------------------------------------------
## there is a bug here if times are integer and constant (or something)
## I think it has to do with boundary.lev creation, as subsequent trips are out of whack

## this fails (but ok if tms is + 1:10)
## d <- data.frame(x=1:10, y=rnorm(10), tms=Sys.time() + c(1:5, 1:5), id=gl(2, 5))
## coordinates(d) <- ~x+y
## tr <- trip(d, c("tms", "id"))

## bound.dates <- seq(min(tr$tms)-1, max(tr$tms)+1, length=5)
## trip.list <- trip.split.exact(tr, bound.dates)

##' @importFrom maptools spRbind
##' @rdname trip-internal
.tripRbind <- function (obj, x) {
    ## not needed, and not possible since classes imported using
    ## NAMESPACE MDS 2012-10-09
    ## suppressMessages(require(maptools))
    tor1 <- getTORnames(obj)
    tor2 <- getTORnames(x)
    if (! identical(tor1, tor2)) stop("trips are not equivalent for rbind")
    SP <- spRbind(as(obj, "SpatialPoints"), as(x, "SpatialPoints"))
    df <- rbind(slot(obj, "data"), slot(x, "data"))
    dupes <- duplicated(cbind(coordinates(SP), df))
    x <- SpatialPointsDataFrame(SP, data=df)[!dupes, ]
    trip(x, tor1)
}

##' @rdname trip-internal
.single.trip.split <- function(tr1, boundary.dates) {
    diff.d <- diff(unclass(boundary.dates))
    if (any(diff.d < 0))
        stop("boundary dates must must sort increasingly")
    if (any(!diff.d > 0))
        stop("duplicates in boundary dates")
    tor <- getTORnames(tr1)
    ## single id trip object
    x <- tr1[, tor]
    x <- data.frame(coordinates(x), x@data[,tor])
    if (min(boundary.dates) > min(x[, 3]))
        stop("boundary dates do not encompass trip range (MIN)")
    if (max(boundary.dates) < max(x[, 3]))
        stop("boundary dates do not encompass trip range (MAX)")
    which.dates <- boundary.dates[boundary.dates > min(x[, 3]) &
                                  boundary.dates < max(x[, 3])]
    which.dates <- rep(which.dates, each=2)
    if (!length(which.dates) > 0) {
        ## we are done
        tr1$boundary.lev <- 1
        res <- list(tr1)
        ind <- which.min(boundary.dates < min(x[, 3]) )
        boundary.names <- paste(boundary.dates[c(ind - 1, ind)],
                                collapse=" <-> ")
        names(res) <- boundary.names
        return(res)
    }
    boundary.ids <- which(boundary.dates > min(x[, 3]) &
                          boundary.dates < max(x[, 3]))
    boundary.ids <- c(boundary.ids[1] - 1, boundary.ids,
                      boundary.ids[length(boundary.ids)] + 1)
    boundary.names <- paste(boundary.dates[boundary.ids[-length(boundary.ids)]],
                            boundary.dates[boundary.ids[-1]], sep=" <-> ")
    fx <- approxfun(x[, 3], x[, 1])
    fy <- approxfun(x[, 3], x[, 2])
    new.x <- fx(which.dates)
    new.y <- fy(which.dates)
    new.1 <- data.frame(new.x, new.y, which.dates,
                        rep(x[1, 4], length(which.dates)))
    names(new.1) <- names(x)
    x.new <- rbind(x, new.1)
    ## sort records
    x.new <- x.new[order(x.new[, 3]), ]
    edges <- which(x.new[, 3] %in% which.dates)
    ## boundary.lev
    boundary.lev <- cumsum(x.new[, 3] %in% which.dates)
    boundary.lev[boundary.lev %% 2 > 0] <-
        boundary.lev[boundary.lev %% 2 > 0] - 1
    x.new$boundary.lev <- unclass(factor(boundary.lev))
    t.list <- split(x.new, x.new$boundary.lev)
    if (!length(t.list) == length(boundary.names))
        stop("names and split do not match")
    names(t.list) <- boundary.names
    ## deal with trips that are too short
    for (i in 1:length(t.list)) {
        if (nrow(t.list[[i]]) < 2)
            stop("this should never happen")
        if (nrow(t.list[[i]]) < 3) {
            x <- t.list[[i]]
            fx <- approxfun(x[, 3], x[, 1])
            fy <- approxfun(x[, 3], x[, 2])
            which.dates <- seq(min(x[, 3]), max(x[, 3]), length=3)
            x1 <- data.frame(fx(which.dates), fy(which.dates),
                             which.dates, rep(x[1, 4], 3), rep(x[1, 5], 3))
            names(x1) <- names(x)
            t.list[[i]] <- x1
        }
    }
    res <- lapply(t.list, function(x) {
        SpatialPointsDataFrame(as.matrix(x[, 1:2]),
                               x[, -c(1, 2)],
                               proj4string=CRS(proj4string(tr1)))
    })
                                        #browser()
    lapply(res, trip, tor)
}




#'
#' Split trip events into exact time-based boundaries.
#'
#'
#' Split trip events within a single object into exact time boundaries, adding
#' interpolated coordinates as required.
#'
#'
#' Motion between boundaries is assumed linear and extra coordinates are added
#' at the cut points.
#'
#' @param x A trip object.
#' @param breaks A character string such as the \code{breaks} argument
#' for \code{\link{cut.POSIXt}}, or alternatively a vector of
#' date-time boundaries. (If the latter hese must encompass all the time range of
#' the entire trip object.)
#' @param \dots Unused arguments.
#' @return
#'
#' A list of trip objects, named by the time boundary in which they lie.
#' @author Michael D. Sumner and Sebastian Luque
#' @details This function was completely rewritten in version 1.1-20. 
#' @seealso See also \code{\link{tripGrid}}.
#' @keywords manip chron
#' @examples
#'
#' \dontrun{
#' set.seed(66)
#' d <- data.frame(x=1:100, y=rnorm(100, 1, 10),
#'                 tms= as.POSIXct(as.character(Sys.time()), tz = "GMT") + c(seq(10, 1000, length=50),
#'                 seq(100, 1500, length=50)), id=gl(2, 50))
#' coordinates(d) <- ~x+y
#' tr <- trip(d, c("tms", "id"))
#'
#' cut(tr, "200 sec")
#'
#' bound.dates <- seq(min(tr$tms) - 1, max(tr$tms) + 1, length=5)
#' trip.list <- cut(tr, bound.dates)
#' bb <- bbox(tr)
#' cn <- c(20, 8)
#' g <- GridTopology(bb[, 1], apply(bb, 1, diff) / (cn - 1), cn)
#'
#' tg <- tripGrid(tr, grid=g)
#' tg <- as.image.SpatialGridDataFrame(tg)
#' tg$x <- tg$x - diff(tg$x[1:2]) / 2
#' tg$y <- tg$y - diff(tg$y[1:2]) / 2
#'
#' op <- par(mfcol=c(4, 1))
#' for (i in 1:length(trip.list)) {
#'   plot(coordinates(tr), pch=16, cex=0.7)
#'   title(names(trip.list)[i], cex.main=0.9)
#'   lines(trip.list[[i]])
#'   abline(h=tg$y, v=tg$x, col="grey")
#'   image(tripGrid(trip.list[[i]], grid=g), interpolate=FALSE,
#'   col=c("white", grey(seq(0.2, 0.7,  length=256))),add=TRUE)
#'   abline(h=tg$y, v=tg$x,  col="grey")
#'   lines(trip.list[[i]])
#'   points(trip.list[[i]], pch=16, cex=0.7)
#' }
#'
#' par(op)
#' print("you may need to resize the window to see the grid data")
#'
#' cn <- c(200, 80)
#' g <- GridTopology(bb[, 1], apply(bb, 1, diff) / (cn - 1), cn)
#'
#' tg <- tripGrid(tr, grid=g)
#' tg <- as.image.SpatialGridDataFrame(tg)
#' tg$x <- tg$x - diff(tg$x[1:2]) / 2
#' tg$y <- tg$y - diff(tg$y[1:2]) / 2
#'
#' op <- par(mfcol=c(4, 1))
#' for (i in 1:length(trip.list)) {
#'   plot(coordinates(tr), pch=16, cex=0.7)
#'   title(names(trip.list)[i], cex.main=0.9)
#'   image(tripGrid(trip.list[[i]], grid=g, method="density", sigma=1),
#'         interpolate=FALSE,
#'         col=c("white", grey(seq(0.2, 0.7, length=256))),
#'         add=TRUE)
#'   lines(trip.list[[i]])
#'   points(trip.list[[i]], pch=16, cex=0.7)
#' }
#'
#' par(op)
#' print("you may need to resize the window to see the grid data")
#'
#' }
#'
#' @method cut trip
#' @export
cut.trip <-
function (x, breaks, ...)
{
    if ("dates" %in% names(list(...))) warning("please use \'breaks\' not \'dates\'")
  tor <- getTORnames(x)
    if (is.character(breaks)) {
        if (length(breaks) > 1) stop("if breaks is character, length(breaks) should be 1L")
        levs <- levels(cut(x[[tor[1]]], breaks))
        datebounds <- seq(as.POSIXct(levs[1L], tz = "GMT"), by = breaks, length = length(levs) + 1)
        breaks <- datebounds
    }
    
   
      uid <- unique(x[[tor[2]]])
      l <- vector("list", length(uid))
      for (i in seq_along(l)) l[[i]] <- cut.one.trip(x[x[[tor[2]]] == uid[i], ], breaks)
      l2 <- vector("list", length(l[[1]]))
      for (j in seq_along(l2)) l2[[j]] <- do.call(rbind, lapply(l, function(x) x[[j]]) )
      lapply(l2[!sapply(l2, is.null)], function(xx) trip(SpatialPointsDataFrame(SpatialPoints(as.matrix(xx[,1:2]), proj4string = CRS(proj4string(x))), xx[,-c(1, 2)]), c("time", "id")))
      
    }
    
    
    

    cut.one.trip <- function(x, breaks, ...) {
      tor <- getTORnames(x)
      coords <- coordinates(x)
      time <- x[[tor[1]]]
      id <- x[[tor[2]]]
      unbrks <- as.numeric(breaks)
      untime <- as.numeric(time)
      no_insert <- unbrks %in% untime
      unbrks2 <- unbrks[!no_insert]
      newbreaks <- sort(c(unbrks[-c(1, length(unbrks))], unbrks2[-c(1, length(unbrks2))], untime))
      newx <- approxfun(untime, coords[,1], rule = 2)(newbreaks)
      newy <- approxfun(untime, coords[,2], rule = 2)(newbreaks)
      
      ntrack <- list(x = newx, y = newy, time = ISOdatetime(1970, 1, 1, 0, 0, 0, tz = "UTC") + newbreaks, 
                     id = rep(id[1], length = length(newx)))
      
      out <- split(as.data.frame(ntrack), cumsum(duplicated(newbreaks)))
      short <- sapply(out, nrow) < 3
      out <- lapply(out, function(x) if (nrow(x) < 3) NULL else x)
      
      for (i in seq_along(out)[-1]) {
        if (short[i] & !short[i-1]) last <- tail(out[[i-1]], 1)
        if (!short[i] & short[i-1] & exists("last")) out[[i]][1,] <- last
      }
      out
    }
    
    

    