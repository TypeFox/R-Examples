#' Filter track data for speed
#' 
#' 
#' Create a filter of a track for "bad" points implying a speed of motion that
#' is unrealistic.
#' 
#' 
#' Using an algorithm (McConnnell et al, 1992), points are tested for speed
#' between previous / next and 2nd previous / next points.  Contiguous sections
#' with an root mean square speed above a given maximum have their highest rms
#' point removed, then rms is recalculated, until all points are below the
#' maximum.  By default an (internal) root mean square function is used, this
#' can be specified by the user.
#' 
#' If the coordinates of the \code{trip} data are not projected, or NA the
#' distance calculation assumeds longlat and kilometres (great circle). For
#' projected coordinates the speed must match the units of the coordinate
#' system.  (The PROJ.4 argument "units=km" is suggested).
#' 
#' @param x trip object
#' @param max.speed speed in kilometres per hour
#' @param test cut the algorithm short and just return first pass
#' @return
#' 
#' Logical vector matching positions in the coordinate records that pass the
#' filter.
#' @note
#' 
#' This algorithm was originally taken from IDL code by David Watts at the
#' Australian Antarctic Division, and used in various other environments before
#' the development of this version.
#' @section Warning:
#' 
#' This algorithm is not considered to be particularly relevant to the problems
#' involved with location uncertainty in animal tracking.  It is provided
#' merely as an illustrative benchmark for further work.
#' 
#' It is possible for the filter to become stuck in an infinite loop, depending
#' on the function passed to the filter.  Several minutes is probably too long
#' for hundreds of points, test on smaller sections if unsure.
#' @author David Watts and Michael D. Sumner
#' @seealso \code{\link{trip}}
#' @references
#' 
#' The algorithm comes from McConnell, B. J. and Chambers, C. and Fedak, M. A.
#' (1992) Foraging ecology of southern elephant seals in relation to the
#' bathymetry and productivity of the southern ocean.  Antarctic Science
#' \emph{4} 393-398
#' @keywords manip
#' @export speedfilter
speedfilter <- function (x, max.speed=NULL, test=FALSE) {
    if (!is(x, "trip"))
        stop("only trip objects supported")
    projected <- is.projected(x)
    if (is.na(projected)) {
        projected <- FALSE
        warning("coordinate system is NA, assuming longlat . . .")
    }
    if (is.null(max.speed)) {
        print("no max.speed given, nothing to do here")
        return(x)
    }
    longlat <- !projected
    coords <- coordinates(x)
    tids <- getTimeID(x)
    time <- tids[, 1]
    id <- factor(tids[, 2])
    x <- coords[, 1]
    y <- coords[, 2]
    pprm <- 3
    grps <- levels(id)
    if (length(x) != length(y))
        stop("x and y vectors must be of same\nlength")
    if (length(x) != length(time))
        stop("Length of times not equal to number of points")
    okFULL <- rep(TRUE, nrow(coords))
    if (test)
        res <- list(speed=numeric(0), rms=numeric(0))
    for (sub in grps) {
        ind <- id == sub
        xy <- matrix(c(x[ind], y[ind]), ncol=2)
        tms <- time[ind]
        npts <- nrow(xy)
        if (pprm%%2 == 0 || pprm < 3) {
            msg <- paste("Points per running mean should be odd and",
                         "greater than 3, pprm=3")
            stop(msg)
        }
        RMS <- rep(max.speed + 1, npts)
        offset <- pprm - 1
        ok <- rep(TRUE, npts)
        if (npts < (pprm + 1)) {
            warning("Not enough points to filter ID: \"", sub,
                    "\"\n continuing . . . \n")
            okFULL[ind] <- ok
            next
        }
        index <- 1:npts
        iter <- 1
        while (any(RMS > max.speed, na.rm=TRUE)) {
            n <- length(which(ok))
            x1 <- xy[ok, ]
            speed1 <- trackDistance(x1[-nrow(x1), 1], x1[-nrow(x1), 2],
                                    x1[-1, 1], x1[-1, 2],
                                    longlat=!projected) /
                                        (diff(unclass(tms[ok])) / 3600)
            speed2 <- trackDistance(x1[-((nrow(x1) - 1):nrow(x1)), 1],
                                    x1[-((nrow(x1) - 1):nrow(x1)), 2],
                                    x1[-(1:2), 1], x1[-(1:2), 2],
                                    longlat=!projected) /
                                        ((unclass(tms[ok][-c(1, 2)]) -
                                          unclass(tms[ok][-c(n - 1, n)])) /
                                         3600)
            thisIndex <- index[ok]
            npts <- length(speed1)
            if (npts < pprm)
                next
            sub1 <- rep(1:2, npts - offset) + rep(1:(npts - offset), each=2)
            sub2 <- rep(c(0, 2), npts - offset) +
                rep(1:(npts - offset), each=2)
            rmsRows <- cbind(matrix(speed1[sub1], ncol=offset, byrow=TRUE),
                             matrix(speed2[sub2], ncol=offset, byrow=TRUE))
            RMS <- c(rep(0, offset),
                     sqrt(rowSums(rmsRows ^ 2) / ncol(rmsRows)))
            if (test & iter == 1) {
                res$speed <- c(res$speed, 0, speed1)
                res$rms <- c(res$rms, 0, RMS)
                break
            }
            RMS[length(RMS)] <- 0
            bad <- RMS > max.speed
            segs <- cumsum(c(0, abs(diff(bad))))
            ## try wrapping ifelse here? no index is quicker
            rmsFlag <- unlist(lapply(split(RMS, segs), function(x) {
                ifelse((1:length(x)) == which.max(x), TRUE, FALSE)
            }), use.names=FALSE)
            rmsFlag[!bad] <- FALSE
            RMS[rmsFlag] <- -10
            ok[thisIndex][rmsFlag > 0] <- FALSE
        }
        okFULL[ind] <- ok
    }
    if (test)
        return(res)
    okFULL
}


# $Id: filter.penSS.R 68 2013-03-20 03:11:06Z sluque $



#' Non-destructive smoothing filter
#' 
#' 
#' Non-destructive filter for track data using penalty smoothing on velocity.
#' 
#' 
#' Destructive filters such as \code{\link{speedfilter}} can be recast using a
#' penalty smoothing approach in the style of Green and Silverman (1994).
#' 
#' This filter works by penalizing the fit of the smoothed track to the
#' observed locations by the sum of squared velocities.  That is, we trade off
#' goodness of fit against increasing the total sum of squared velocities.
#' 
#' When lambda=0 the smoothed track reproduces the raw track exactly.
#' Increasing lambda favours tracks requiring less extreme velocities, at the
#' expense of reproducing the original locations.
#' @name filter.penSS
#' @param tr A \code{trip} object.
#' @param lambda Smoothing parameter, see Details.
#' @param first Fix the first location and prevent it from being updated by the
#' filter.
#' @param last Fix the last location and prevent it from being updated by the
#' filter.
#' @param \dots Arguments passed on to \code{\link{nlm}}
#' @return
#' 
#' A trip object with updated coordinate values based on the filter - all the
#' data, including original coordinates which are maintained in the trip data
#' frame.
#' @author Simon Wotherspoon and Michael Sumner
#' @seealso \code{\link{speedfilter}}
#' @references
#' 
#' Green, P. J. and Silverman, B. W. (1994). Nonparametric regression and
#' generalized linear models: a roughness penalty approach. CRC Press.
#' @keywords manip misc
#' @examples
#' 
#' 
#' \dontrun{## Example takes a few minutes
#' 
#' ## Fake some data
#' 
#' ## Brownian motion tethered at each end
#' brownian.bridge <- function(n, r) {
#'   x <- cumsum(rnorm(n, 0, 1))
#'   x <- x - (x[1] + seq(0, 1, length=n) * (x[n] - x[1]))
#'   r * x
#' }
#' 
#' ## Number of days and number of obs
#' days <- 50
#' n <- 200
#' 
#' ## Make separation between obs gamma distributed
#' x <- rgamma(n, 3)
#' x <- cumsum(x)
#' x <- x/x[n]
#' 
#' ## Track is lissajous + brownian bridge
#' b.scale <- 0.6
#' r.scale <- sample(c(0.1, 2, 10.2), n, replace=TRUE,
#'                   prob=c(0.8, 0.18, 0.02))
#' set.seed(44)
#' 
#' tms <- ISOdate(2001, 1, 1) + trunc(days * 24 * 60 * 60 *x)
#' lon <- 120 + 20 * sin(2 * pi * x) +
#'     brownian.bridge(n, b.scale) + rnorm(n, 0, r.scale)
#' lat <- -40 + 10 *(sin(3 * 2 * pi * x) + cos(2 * pi * x) - 1) +
#'     brownian.bridge(n, b.scale) + rnorm(n, 0, r.scale)
#' 
#' tr <- new("trip",
#'           SpatialPointsDataFrame(cbind(lon, lat),
#'                                  data.frame(gmt=tms, id="lbb")),
#'                                  TimeOrderedRecords(c("gmt", "id")))
#' plot(tr)
#' 
#' ## the filtered version
#' trf <- filter.penSS(tr, lambda=1, iterlim=400, print.level=1)
#' 
#' lines(trf)
#' 
#' }
#' 
#' 
#' @export filter.penSS
filter.penSS <- function(tr, lambda, first=TRUE, last=TRUE,...) {

    penalized <- function(x) {
        ## Form smoothed track
        p <- p.obs
        p[sub, ] <- x
        ## Velocities between smoothed points
        ##v <- gc.dist(p[-n,],p[-1,])/dt
        v <- trackDistance(p[, 2:1]) / dt
        ## Distances from smoothed points to observations
        ##d <- gc.dist(p,p.obs)
        d <- trackDistance(p[, 2], p[, 1], p.obs[, 2], p.obs[, 1])
        ## This is the penalized sum of squares
        (sum(d ^ 2) + lambda * sum(v ^ 2)) / n ^ 2
    }

    if (length(summary(tr)$tripID) > 1) {
        msg <- paste("trip object contains multiple events,",
                     "only the first trip used")
        warning(msg)
        tr <- tr[tr[[getTORnames(tr)[2]]] == summary(tr)$tripID[1], ]
    }
    ## Number of points and subset
    n <- nrow(tr)
    sub <- (1 + first):(n - last)
    ## Observed points
    ##  p.obs <- as.matrix(tr[,c("Lat","Lon")])
    p.obs <- coordinates(tr)[, 2:1]
    ## Time intervals (in days) between obs
    ##dt <- diff(unclass(tr$Time)/(24*60*60))
    dt <- diff(unclass(tr[[getTORnames(tr)[1]]]) / (24 * 60 * 60))
    mn <- nlm(penalized, as.matrix(p.obs[sub, ]), ...)
    m <- n - (first + last)
    res <- coordinates(tr)
    ##  tr$Lat[sub] <- mn$estimate[1:m]
    ##  tr$Lon[sub] <- mn$estimate[m+1:m]
    res[sub, 2] <- mn$estimate[1:m]
    res[sub, 1] <- mn$estimate[m + 1:m]
    res <- SpatialPointsDataFrame(res, as.data.frame(tr),
                                  proj4string=CRS(proj4string(tr)))
    trip(res, getTORnames(tr))
}



###_ + Emacs local variables
## Local variables:
## allout-layout: (+ : 0)
## End:
