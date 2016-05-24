#############################################################################
## package 'secr'
## methods.R
## Methods for classes traps, capthist and mask
## Last changed
## 2010 04 27 version 1.4
## 2010 04 30 scrapped read.captures
## 2010 05 02 usage may be separated by commas in read.traps
## 2010 05 03 bugfix in subset.capthist for 'multi' data
## 2010 06 08 sim.popn moved to own file sim.popn.R
## 2010 07 01 alphanumeric detectfn option
## 2010 06 08 fixed bug in insidepoly (wrong call of C fn)
## 2010 07 03 rbind.mask handles covariates
## 2010 08 30 print.mask removed
## 2010 10 10 enabled compound halfnormal in secr.fit
## 2010 10 10 added function 'detectpar'
## 2010 10 11 secr.fit default binomN changed to Poisson if detector == count
## 2010 10 12 covariates<- modified for multi-session data
## 2010 10 12 plot.mask now handles multi-session mask
## 2010 10 15 function ms to recognise multisession objects
## 2010 10 15 make.mask poly clipping extended to all types
## 2010 10 19 version 1.5
## 2010-10-21 AIC.secrlist
## 2010-10-22 predict.secrlistyes## 2010-10-24 multisession detectpar
## 2010-10-25 revision of all error and warning messages in R code
## 2010-11-04 logLik.secr
## 2011-01-04 G&R parameterisation
## 2011-01-12 cue detector  -- removed 2015-10-01
## 2011-02-01 remove binomN from traps object attribute list
## 2011-02-07 major revisions for polygonX, transectX
## 2011-03-18 unmarked detector
##            for counts of unidentified animals at points
## 2011-06-07 rbind.capthist allows for xy, signal
## 2011-06-17 moved insidepoly to utility.R (renamed to pointsInPolygon 2011-10-21)
## 2011-06-23 polygon and transect names retain ordering of input dataframe
##            in make.poly and make.transect
## 2011-06-24 fixed transectlength, searcharea
## 2011-07-08 adjust 'class' to 'inherits' in make.mask
## 2011-08-18 check for null rownames in animalID
## 2011-09-11 subset.popn commissioned
## 2011-10-10 make.mask moved to make.mask.R
## 2011-10-10 make.grid, make.poly etc. moved to make.grid.R
## 2011-10-20 predict.secr modified for user-supplied D function
## 2011-10-20 subset.popn gets poly argument
## 2012-05-14 predict.secr patch for missing 'session' in newdata (cf Deb's problem 7/2/12)
## 2012-07-28 subset.capthist substantial revision for signals
## 2012-08-07 spacing methods for traps, mask gets new argument 'recalculate'
## 2012-08-07 spacing attribute updated after subset, split or rbind of traps
## 2012-09-14 split.traps now in its own file
## 2012-09-14 signalframe extract/replace function
## 2012-10-24 predict.secr drops unused columns of newdata from names of output components
## 2012-10-25 plot.capthist checks for zero captures and proximity tracks
## 2012-12-22 extended usage<-
## 2013-01-16 fixed rbind.traps for usage
## 2013-02-07 timevaryingcov<- bug fixed
## 2013-03-14 predict.secr acquires argument 'type'
## 2013-06-08 print.secr includes Mixture (hcov) :
## 2013-06-10 plot.capthist has safe on.exit return to old palette
## 2013-06-15 shift.mask
## 2013-08-17 subset.capthist hardened so copes with signal detectors with zero rows
## 2013-08-30 detectpar becomes method
## 2013-10-29 fixpmix code for predict.secr moved to utility.R
## 2013-11-20 shifted plot.capthist to separate file
## 2013-12-15 maskarea function for nrow(mask) * attr(mask,'area')
## 2014-02-08 trim.secr() added call to default vector of components dropped
## 2014-02-19 secrlist '[' method
## 2014-05-31 summary.capthist counts losses when dim>2
## 2014-05-31 read.mask spurious error msg with columns argument
## 2014-08-19 predict.secr transferred to predict.secr.r
## 2014-09-06 summary.mask slightly rearranged, and allows linearmask
## 2014-09-14 plot.mask moved to plot.mask.r
## 2014-12-24 print.summary.capthist allows for nonspatial capthist
## 2015-01-16 print.secr shows userdist type
## 2015-03-31 parameter count allows for details$miscparm
## 2015-09-02 change alongtransect to fix occsim: remove 'S' prefix
## 2015-10-02 summary.capthist reports sightings
## 2015-10-10 subset.capthist updated for mark-resight data
###############################################################################

# Generic methods for extracting attributes etc

usage      <- function (object, ...) UseMethod("usage")
markocc    <- function (object, ...) UseMethod("markocc")
Tu         <- function (object, ...) UseMethod("Tu")
Tm         <- function (object, ...) UseMethod("Tm")
clusterID  <- function (object, ...) UseMethod("clusterID")
clustertrap <- function (object, ...) UseMethod("clustertrap")
covariates <- function (object, ...) UseMethod("covariates")
traps      <- function (object, ...) UseMethod("traps")
detector   <- function (object, ...) UseMethod("detector")
spacing    <- function (object, ...) UseMethod("spacing")
session    <- function (object, ...) UseMethod("session")
trim       <- function (object, drop, keep) UseMethod("trim")
timevaryingcov <- function (object, ...) UseMethod("timevaryingcov")

rotate     <- function (object, degrees, centrexy=NULL, ...) UseMethod("rotate")
shift      <- function (object, shiftxy, ...) UseMethod("shift")
flip       <- function (object, lr=F, tb=F, ...) UseMethod("flip")

ms         <- function (object, ...) UseMethod("ms")
detectpar  <- function (object, ...) UseMethod("detectpar")
signal     <- function (object, ...) UseMethod("signal")
noise      <- function (object, ...) UseMethod("noise")

# Default methods for specialised functions

ms.default <- function (object, ...)       {
    inherits(object, 'list')
}
ms.mask <- function (object, ...)       {
    !is.data.frame(object)
}
ms.secr <- function (object, ...)       {
    ms(object$capthist)
}

usage.default <- function (object, ...)       {
    if (ms(object)) lapply(object, usage.default, ...)
    else attr(object,'usage')
}

markocc.default <- function (object, ...)       {
    if (ms(object)) lapply(object, markocc.default, ...)
    else attr(object,'markocc')
}

Tu.default <- function (object, ...)       {
    if (ms(object)) lapply(object, Tu.default, ...)
    else attr(object,'Tu')
}

Tm.default <- function (object, ...)       {
    if (ms(object)) lapply(object, Tm.default, ...)
    else attr(object,'Tm')
}

sighting <- function(object) {
    markocc <- markocc(object)
    if (is.null(markocc)) 
        FALSE
    else 
        any(markocc < 1) 
}
    

clusterID.default <- function (object, ...)       {
    if (ms(object)) lapply(object, clusterID.default, ...)
    else attr(object,'cluster')
}

clustertrap.default <- function (object, ...)       {
    if (ms(object)) lapply(object, clustertrap.default, ...)
    else attr(object,'clustertrap')
}

covariates.default <- function (object, ...)  {
    if (ms(object)) lapply(object, covariates.default, ...)
    else attr(object,'covariates')
}

timevaryingcov.default <- function (object, ...)  {
    if (ms(object)) lapply(object, timevaryingcov.default, ...)
    else attr(object,'timevaryingcov')
}

traps.default <- function (object, ...)       {
    if (ms(object)) {
        temp <- lapply(object, traps.default, ...)
        class(temp) <- c('list', 'traps')
        temp
    }
    else{
        attr(object,'traps')
    }
}

detector.default <- function (object, ...)    {
## assumed constant across MS
    if (ms(object)) {
        detector.default (object[[1]], ...)
    }
    else {
        if (is.null(object)) NULL
        else attr(object,'detector')
    }
}

spacing.default <- function (object, ...)    {
    if (is.null(object))
        NULL
    else {
        attr(object,'spacing')
    }
}

spacing.traps <- function (object, ..., recalculate = FALSE)    {
    if (ms(object)) {
        sapply(object, spacing.traps, ..., recalculate)
    }
    else {
        if (is.null(object)) {
            NULL
        }
        else {
            temp <- attr(object,'spacing')
            if ((is.null(temp) | recalculate) & (nrow(object)>1)) {
                spacing <- as.matrix(dist(object))
                sp <- apply(spacing,1,function(x) min(x[x>0]))
                mean(sp)
            }
            else
                temp
        }
    }
}

spacing.mask <- function (object, ..., recalculate = FALSE)    {
    if (ms(object)) {
        sapply(object, spacing.mask, ...)
    }
    else {
        if (is.null(object)) NULL
        else {
            temp <- attr(object,'spacing')
            if ((is.null(temp) | recalculate)& (nrow(object)>1) ) {
                spacing <- as.matrix(dist(object))
                sp <- apply(spacing,1,function(x) min(x[x>0]))
                mean(sp)
            }
            else
                temp
        }
    }
}

## traps object
polyID <- function (object)    {
    if (ms(object)) {
        polyID (object[[1]])
    }
    else {
        if (inherits(object,'traps')) {
            temp <- attr(object,'polyID')
#            if (is.null(temp)) temp <- factor(rep(1,nrow(object)))
            if (is.null(temp)) temp <- factor(1:nrow(object))   ## all different
            temp
        }
        else
        if (inherits(object,'capthist')) {
            stop ("use trap() to extract polyID from 'capthist' object")
        }
        else stop ("polyID requires 'traps' object")
    }
}

## traps object
transectID <- function (object)    {
    if (ms(object)) {
        transectID (object[[1]])
    }
    else {
        if (inherits(object,'traps')) {
            if (!detector(object) %in% c('transect','transectX'))
                stop ("requires transect detector")
            temp <- attr(object,'polyID')
#            if (is.null(temp)) temp <- factor(rep(1,nrow(object)))
            if (is.null(temp)) temp <- factor(1:nrow(object))
            temp
        }
        else
        if (inherits(object,'capthist')) {
            stop ("use trap() to extract transectID from 'capthist' object")
        }
        else stop ("transectID requires 'traps' object")
    }
}

xy <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")

    if (ms(object)) {
        lapply(object, xy)
    }
    else {
        if (detector(traps(object)) %in%
            c('polygonX', 'transectX', 'polygon','transect','telemetry')) {
            attr(object, 'detectedXY')
        }
        else
            NULL
    }
}

telemetryxy <- function (object, includeNULL = FALSE) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")

    if (ms(object)) {
        lapply(object, telemetryxy, includeNULL)
    }
    else {
        xylist <- attr(object, 'xylist')
        if (includeNULL) {
            ## expand for untelemetered animals in object
            nullxy <- rep(list(matrix(nrow=0, ncol=2)), nrow(object))
            names(nullxy) <- rownames(object)
            xylist <- c(xylist, nullxy[!(rownames(object) %in% names(xylist))])
        }
        xylist
    }
}

telemetered <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    if (ms(object)) {
        lapply(object, telemetered)
    }
    else {
        rownames(object) %in% names(telemetryxy(object))
    }
}

alongtransect <- function (object, tol = 0.01) {
    ptalongtransect <- function (i) {
        ## where is point i on its transect k?

        k <- trans[i]
        transectxy <- as.matrix(lxy[[k]])
        nr <- nrow(transectxy)
        temp <- .C('alongtransect',  PACKAGE = 'secr',
            as.double (xyi[i,]),
            as.integer (0),
            as.integer (nr-1),
            as.integer (nr),
            as.double (transectxy),
            as.double (tol),
            result = double(1))
        temp$result
    }
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")

    if (ms(object)) {
        lapply(object, alongtransect)
    }
    else {
        trps <- traps(object)
        
        if (detector(trps) %in% c('transectX', 'transect')) {
            
#             trans <- trap(object, names = TRUE)
#             xyi <- xy(object)
#             ## 2015-09-02 change to fix occsim: remove 'S' prefix
#             lxy <- split (trps, levels(transectID(trps)), prefix = "")
            
            trans <- trap(object, names = FALSE)
            xyi <- xy(object)
            lxy <- split (trps, transectID(trps))

            sapply(1:nrow(xyi), ptalongtransect)
        }
        else
            NULL
    }
}

clusterID <- function (object) {
    if (ms(object)) {
        lapply(object, clusterID)
    }
    else {
        if (inherits(object, 'capthist')) {
            trps <- traps(object)
            clusterID(trps)[trap(object, names=FALSE)]
        }
        else
            attr(object, 'cluster')
    }
}

clustertrap <- function (object) {
    if (ms(object)) {
        lapply(object, clustertrap)
    }
    else {
        if (inherits(object, 'capthist')) {
            trps <- traps(object)
            clustertrap(trps)[trap(object, names=FALSE)]
        }
        else
        attr(object, 'clustertrap')
    }
}

signal.default <- function(object, ...) {
    stop ("only for capthist data")
}

signal.capthist <- function (object, ...) {
    if (ms(object)) {
        lapply(object, signal)
    }
    else {
        if (detector(traps(object)) %in% c('signal','signalnoise')) {
            attr(object, 'signalframe')$signal
        }
        else
            NULL
    }
}

noise.default <- function(object, ...) {
    stop ("only for capthist data")
}

noise.capthist <- function (object, ...) {
    if (ms(object)) {
        lapply(object, noise)
    }
    else {
        if (detector(traps(object)) %in% c('signalnoise')) {
            attr(object, 'signalframe')$noise
        }
        else
            NULL
    }
}

signalframe <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")

    if (ms(object)) {
        lapply(object, signalframe)
    }
    else {
        if (detector(traps(object)) %in% c('signal','signalnoise')) {
            attr(object, 'signalframe')
        }
        else
            NULL
    }
}

signalmatrix <- function (object, noise = FALSE, recodezero = FALSE,
    prefix = 'Ch', signalcovariates = NULL, names = NULL) {
    if (ms(object)){
        lapply(object, signalmatrix, noise = noise, recodezero = recodezero,
               prefix=prefix, signalcovariates = signalcovariates, names = names)
    }
    else {
        if (!detector(traps(object)) %in% .localstuff$detectors3D)
            stop ("requires 3-D detector")

#        reordered<- function(x) {
#           x[order(trap(object), occasion(object), animalID(object))]
#        }
        tmpsignal <- object
        if (dim(tmpsignal)[2]>1)
            warning("using only first occasion")
        if (noise)
            sound <- noise(object)
        else
            sound <- signal(object)
        if (!is.null(sound))
## error found 2014-01-24
#            tmpsignal[(tmpsignal>0) | is.na(tmpsignal)] <- reordered(sound)
            tmpsignal[(tmpsignal>0) | is.na(tmpsignal)] <- sound
        if (recodezero)
            tmpsignal[tmpsignal==0] <- NA
        tmpsignal <- tmpsignal[,1,,drop=FALSE]
        tmpsignal <- as.data.frame(tmpsignal)
        ## added 2013-09-03
        if (!is.null(names)) {
            if (length(names)==ncol(tmpsignal))
                names(tmpsignal) <- names
            else
                names(tmpsignal) <- paste(prefix, 1:ncol(tmpsignal), sep='')
        }
        else {
            names(tmpsignal) <- paste(prefix, rownames(traps(object)), sep='')
        }
        if (!is.null(signalcovariates)) {
            sf <- signalframe(object)
            firstID <- match(rownames(tmpsignal), sf$Selection)
            tmpsignal[,signalcovariates] <- sf[firstID,signalcovariates]
        }
        tmpsignal
    }
}

occasion <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    if (ms(object)) {
        lapply(object, occasion)
    }
    else {
        if (detector(traps(object)) %in% .localstuff$exclusivedetectors) {
            col(object)[abs(object)>0]
        }
        else {
            apo <- aperm(object,c(2,1,3))
            temp <- matrix(apo, nrow = dim(apo)[1])
            temp <- array(row(temp), dim = dim(apo))
            ## rep(aperm(temp,c(2,1,3)), abs(object)) 
            ## 2015-11-03 avoid array when empty
            as.numeric(rep(aperm(temp,c(2,1,3)), abs(object)))
        }
    }
}

alive <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    if (ms(object)) {
        lapply(object, alive)
    }
    else
        if (detector(traps(object)) %in% .localstuff$exclusivedetectors) {
            sign(object)[abs(object)>0] > 0
        }
        else {
            temp <- sign(object)[abs(object)>0] > 0
            rep(temp, abs(object[abs(object)>0]))
        }
}

trap <- function (object, names = TRUE) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    if (ms(object)) {
        lapply(object, trap, names=names)
    }
    else {
        if (names)
            values <- row.names(traps(object))
        else
            values <- 1:nrow(traps(object))
        if (detector(traps(object)) %in% .localstuff$exclusivedetectors) {
            values[abs(object[abs(object)>0])]
        }
        else {
            apo <- aperm(object,c(3,1,2))
            temp <- matrix(apo, nrow = dim(apo)[1])
            temp <- array(row(temp), dim = dim(apo))
            k <- aperm(temp,c(2,3,1))
            ## rep(values[k], abs(object))
            ## 2015-11-03
            if (names) 
                as.character(rep(values[k], abs(object)))
            else
                as.numeric(rep(values[k], abs(object)))
        }
    }
}

animalID <- function (object, names = TRUE) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    if (ms(object)) {
        lapply(object, animalID, names=names)
    }
    else {
        if (nrow(object) == 0) {
           ## ''
            if (names)
           character(0)  ## 2011-04-08
            else
                numeric(0) ## 2013-09-09
       }
        else {
            if (names & !is.null(row.names(object)))  ## 2011-08-18 null check
                values <- row.names(object)
            else
                values <- 1:nrow(object)
            if (detector(traps(object)) %in% .localstuff$exclusivedetectors) {
                detrow <- row(object)[abs(object)>0]
                values[detrow]
            }
            else {
                temp <- matrix(object, nrow = dim(object)[1])
                n <- array(row(temp), dim=dim(object))
                rep(values[n], abs(object))
            }
        }
    }
}

detectionindex <- function (object) {
## detectionindex is non-exported function 2012-02-11
## to which original cell in dim3 capthist object does a detection relate?
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    if (ms(object)) {
        lapply(object, detectionindex)
    }
    else {
        if (nrow(object) == 0)
           numeric(0)
        else {
            if (detector(traps(object)) %in% .localstuff$exclusivedetectors) {
                stop("detectionindex is not intended  for exclusive detectors")
            }
            else {
                x <- object
                x[] <- 1:length(x)
                x[abs(object)==0] <- 0
                x <- aperm(x, c(3,2,1))
                x <- x[x>0]
                count <- aperm(object, c(3,2,1))
                count <- count[abs(count)>0]
                rep(x,count)
            }
        }
    }
}

polyarea <- function (xy, ha = TRUE) {
    if (inherits(xy, 'SpatialPolygons')) {
        if (!requireNamespace('rgeos', quietly = TRUE))
            stop ("package rgeos is required for area of SpatialPolygons")
        temparea <- rgeos::gArea(xy)
    }
    else {
        nr <- length(xy$x)
        ## beat problem with integer overflow 2010-11-21
        xy$x <- as.double(xy$x)
        xy$y <- as.double(xy$y)
        sc <- sum (xy$x[-nr] * xy$y[-1] - xy$x[-1] * xy$y[-nr])
        temparea <- abs(sc/2)
    }
    ifelse (ha, temparea/10000, temparea)
}

searcharea <- function (object)    {
## discontinued use of searchcell attribute 2015-10-19
## uses spacex.spacey, or spacing^2, or 1
## requires traps object
    if (ms(object)) {
        ## 2011-06-24
        lapply(object, searcharea)
    }
    else {
        if (detector(object) %in% c('polygon','polygonX','telemetry')) {
            ## assume poly closed
            sapply(split(object, levels(polyID(object))), polyarea)
        }
        else {
            spx <- attr(object,'spacex')
            spy <- attr(object,'spacey')
            if (is.null(spx) | is.null(spy))
                sp2 <- attr(object,'spacing')^2
            else
                sp2 <- spx * spy
            ifelse (is.null(sp2), 1, sp2 / 10000)
        }
    }
}

transectlength <- function (object)    {
## requires traps object
    if (ms(object)) {
        ## 2011-06-24
        lapply(object, transectlength)
    }
    else {
        if (detector(object) %in% c('transect','transectX')) {
            calclength <- function (xy) {
                nr <- nrow(xy)    ## number of vertices
                segments <- (xy$x[-nr] - xy$x[-1])^2 + (xy$y[-nr] - xy$y[-1])^2
                sum (segments^0.5)
            }
            sapply(split(object, polyID(object)), calclength)
        }
        else NA
    }
}

session.default <- function (object, ...)     {
## bypass session attribute for multi-session objects 2009 12 22
## use names(object) for lists i.e. multi-session objects

    if (ms(object)) {
       temp <- names(object)
    }
    else {
        temp <- attr(object,'session')
        if (is.null(temp)) temp <- 1    ## added 2010-02-03
    }
    names(temp) <- NULL
    as.character(temp)      ## 2010 02 25
}

trim.default <- function (object, drop, keep)     {
## drop unwanted named components of a list
## conservative resolution of conflicts between drop & keep
    objnames <- names(object)
    indices <- 1:length(object)
    if (missing(drop)) drop <- indices  # all! but wait...
    if (missing(keep)) keep <- 0
    if (is.character(keep)) keep <- match(keep, objnames)
    if (is.character(drop)) drop <- match(drop, objnames)
    drop <- drop[drop %in% indices[! (indices %in% keep)]]
    ## by index, so have to work from end
    for (i in sort(drop, decreasing = T)) object[[i]] <- NULL
    object
}

###############################################################################

rotate.default <- function (object, degrees, centrexy=NULL, ...) {

    rotatefn <- function (xy) {
        # about centre
        x <- xy[1] - centrexy[1]
        y <- xy[2] - centrexy[2]
        x2 <- x * cos(theta) + y * sin(theta) + centrexy[1]
        y2 <- - x * sin(theta) + y * cos(theta) + centrexy[2]
        c(x2,y2)
      }

##    if (ms(object)) lapply(object, rotate.default, degrees=degrees, centrexy=centrexy, ...)
##    else

    object <- as.matrix(object)
    if (dim(object)[2] <2)
        stop ("requires at least 2 columns")
    if (abs(degrees)>0) {
        if (is.null(centrexy)) centrexy <- c(0,0)
        theta <- 2*pi*degrees/360 # convert to radians
        temp <- t(apply (object,1,rotatefn))
        if (!is.null(dimnames(object)[[2]]))
            dimnames(temp)[[2]] <- dimnames(object)[[2]][1:2]
        temp
    } else object[,1:2]
}
###############################################################################

shift.default <- function (object, shiftxy, ...) {
##    if (ms(object)) lapply(object, shift.default, shiftxy=shiftxy, ...)
##    else

    object <- as.matrix(object[,1:2])
    object[,1] <- object[,1] + shiftxy[1]
    object[,2] <- object[,2] + shiftxy[2]
    object
}
###############################################################################

flip.default <- function (object, lr=F, tb=F, ...) {
##    if (ms(object)) lapply(object, flip.default, lr=lr, tb=tb, ...)
##    else
    object <- as.matrix(object[,1:2])
    if (is.logical(lr)) {
        if (lr) object[,1] <- 2 * mean(object[,1]) - object[,1]  ## flip about mean
    } else
        if (is.numeric(lr)) object[,1] <- 2*lr-object[,1]  ## flip about lr

    if (is.logical(tb)) {
        if (tb) object[,2] <- 2 * mean(object[,2]) - object[,2]  ## flip about mean
    } else
        if (is.numeric(tb)) object[,2] <- 2*tb-object[,2]  ## flip about tb

    object
}
###############################################################################

# Generic methods for replacing values

'usage<-' <- function (object, value) {
    
    ## extended 2012-12-22
    ## bug : loses session names when called as usage(traps(ovenCH)) <- c(1,10)
    
    if (ms(object)) {
        ## replicate across sessions as required
        if (!is.list(value)) {
            nonmatrix <- length(dim(value)) != 2
            valuelist <- vector(mode='list', length=length(object))
            if (nonmatrix & (length(value)>2)) {
                if (length(value) != (length(valuelist)+1))
                    stop("invalid value vector - how many occasions?")
                for (i in 1:length(valuelist)) valuelist[[i]] <- value[c(1,i+1)]
            }
            else
                for (i in 1:length(valuelist)) valuelist[[i]] <- value
        }
        else
            valuelist <- value
        temp <- mapply('usage<-', object, valuelist, SIMPLIFY = FALSE)
        class(temp) <- c('list','traps')
        temp
    }
    else {
        if (!is.null(value)) {
            nonmatrix <- length(dim(value)) != 2
            if (nonmatrix) {
                value <- unlist(value)
                if (length(value) == 1) {
                    if (is.null(usage(object)))
                        stop("requires number of occasions")
                    else noccasions <- ncol(usage(object))
                }
                else
                    noccasions <- value[2]
                value <- value[1]
            }
            else {
                noccasions <- NULL
                value <- as.matrix(value)  ## coerce dataframe to matrix
            }
            if (is.null(noccasions)) {
                if (is.null(dim(value)))
                    stop("usage value should be matrix or scalar with noccasions")
                else if (nrow(value) != ndetector(object))
                    stop("mismatch between value and number of detectors")
            }
            else {
                ## assume value is scalar
                value <- matrix(value, nrow = ndetector(object), ncol = noccasions)
            }
            if (!is.null(markocc(object))) {
                if (ncol(value) != length(markocc(object)))
                    stop ("replacement value not compatible with existing markocc attribute")
            }
        }
        structure (object, usage = value)
    }
}

'markocc<-' <- function (object, value) {
    
    if (ms(object)) {
        ## replicate across sessions as required
        if (!is.list(value)) {
            valuelist <- vector(mode='list', length=length(object))
            for (i in 1:length(valuelist)) valuelist[[i]] <- value
        }
        else
            valuelist <- value
        temp <- mapply('markocc<-', object, valuelist, SIMPLIFY = FALSE)
        class(temp) <- c('list','traps')
        temp
    }
    else {
        if (!inherits(object, 'traps'))
            stop("markocc attribute can only be assigned for traps objects")
        if (!is.null(value)) {
            if (!is.null(usage(object))) {
                if (length(value) != ncol(usage(object)))
                    stop ("markocc not compatible with existing usage attribute")
            }
            if (!all(value %in% c(-1,0,1)))
                stop ("markocc values must be -1, 0 or 1")
        }
        structure (object, markocc = value)
    }
}

'Tu<-' <- function (object, value) {
    
    if (ms(object)) {
        if (!is.list(value)) {
            stop("replacement of Tu for multisession object requires a list")
        }
        else {
            temp <- mapply('markocc<-', object, value, SIMPLIFY = FALSE)
            class(temp) <- class(object)
            temp
        }
    }
    else {
        if (!is.null(value)) {
            if (is.null(markocc <- markocc(traps(object)))) 
                stop("cannot assign Tu for capthist in which traps has no markocc attribute")
            if (length(value)>1) {
                if (length(dims <- dim(value)) != 2)
                    stop ("require either total count or traps x occasions matrix")
                if (dims[1] != ndetector(traps(object)))
                    stop ("Tu not compatible with traps attribute")
                if (dims[2] != length(markocc))
                    stop ("Tu not compatible with markocc attribute")
            }
            if (any(value<0))
                stop ("sighting counts cannot be negative")
        }
        structure (object, Tu = value)
    }
}

'Tm<-' <- function (object, value) {
    
    if (ms(object)) {
        if (!is.list(value)) {
            stop("replacement of Tm for multisession object requires a list")
        }
        else {
            temp <- mapply('markocc<-', object, value, SIMPLIFY = FALSE)
            class(temp) <- class(object)
            temp
        }
    }
    else {
        if (!is.null(value)) {
            if (is.null(markocc <- markocc(traps(object)))) 
                stop("cannot assign Tm for capthist in which traps has no markocc attribute")
            if (length(value)>1) {
                if (length(dims <- dim(value)) != 2)
                    stop ("require either total count or traps x occasions matrix")
                if (dims[1] != ndetector(traps(object)))
                    stop ("Tm not compatible with traps attribute")
                if (dims[2] != length(markocc))
                    stop ("Tm not compatible with markocc attribute")
            }
            if (any(value<0))
                stop ("sighting counts cannot be negative")
        }
        structure (object, Tm = value)
    }
}

'clusterID<-' <- function (object, value) structure (object, cluster = value)
'clustertrap<-' <- function (object, value) structure (object, clustertrap = value)

'covariates<-' <- function (object, value) {
## modified for multi-session data 2010-10-15
    if (is.null(value))
        structure (object, covariates = NULL)
    else {
        if (ms(object)) {
            if (length(object) != length(value))
                stop ("mismatch between multisession object and covariates")
            value <- lapply(value, as.data.frame)
        }
        else {
            value <- as.data.frame(value)
            nrequired <- nrow(object)
            if (inherits(object, 'traps'))
                if (detector(object) %in% .localstuff$polydetectors)
                    nrequired <- length(levels(polyID(object)))
            if (nrow(value) != nrequired)
                stop ("length of covariate does not match")
        }
        structure (object, covariates = value)
    }
}

'timevaryingcov<-' <- function (object, value) {
## 2012-10-31, modified 2013-02-07
    if (is.null(value))
        structure (object, timevaryingcov = NULL)
    else {
        if (ms(object)) {
            temp <- object
            if (!(is.list(value[[1]]) & (length(value) == length(object))))
                value <- list(value)
            for (i in 1:length(object))
                timevaryingcov(temp[[i]]) <- value[[i]]
            temp
        }
        else {
            if (!inherits(object, 'traps'))
                stop("timevaryingcov is for traps objects")
            if (!is.list(value) | is.null(names(value)))
                stop("value should be a list of one or more named vectors")
            if (!is.null(usage(object))) {
                OK <- sapply(value, function(x)
# bug fixed 2013-02-07       (ncol(usage) == length(x)))
                             (ncol(usage(object)) == length(x)))
                if (any(!OK))
                    warning ("mismatch between number of occasions in usage and timevaryingcov")
            }
            if (is.character(value))
                if (!all(value %in% names(covariates(object))))
                    warning ("mismatch between character vector and covariate names")
            structure (object, timevaryingcov = value)
        }
    }
}

'detector<-' <- function (object, value) {
    if (!(value %in% .localstuff$validdetectors))
        stop ("invalid detector type")
    structure (object, detector = value)
}

'spacing<-' <- function (object, value) {
    if (!(is.numeric(value)))
        stop ("non-numeric spacing")
    if (ms(object)) {
        stop ("not sure how to replace spacing of ms object")
    }
    else {
        structure (object, spacing = value)
    }
}

'polyID<-' <- function (object, value) {
    if (!inherits(object,'traps'))
        warning ("polyID requires 'traps' object")
    if (length(value)==1) value <- rep(value, nrow(object))
    value <- factor(value)
    structure (object, polyID = value)
}

'transectID<-' <- function (object, value) {
    if (length(value)==1) value <- rep(value, nrow(object))
    value <- factor(value)
    structure (object, polyID = value)
}

'xy<-' <- function (object, value) {
    if (!is.null(value)) {
        if (detector(traps(object)) %in% .localstuff$exclusivedetectors)
            ndetections <- sum(abs(object)>0)
        else
            ndetections <- sum(abs(object))
        if (nrow(value) != ndetections)
            stop ("requires one location per detection")
        if (!(detector(traps(object)) %in% .localstuff$polydetectors) |
                !(inherits(object,'capthist')))
            stop ("requires 'capthist' object with ",
                  "polygon-like detector")
        if (ms(object))
            stop ("requires single-session 'capthist' object")
    }
    structure (object, detectedXY = value)
}

'clusterID<-' <- function (object, value) {

    if (ms(object))
        stop ("cluster requires single-session 'traps' object")

    if (length(value)==1) value <- rep(value, nrow(object))

    if (!(inherits(object, 'traps')))
        stop ("requires clustered 'traps' object")

    if (length(value) > 0) {
        if (length(value) != nrow(object))
            stop ("requires one cluster per detector or detector vertex")
    }
    else
        value <- NULL

    structure (object, cluster = factor(value))
}

'clustertrap<-' <- function (object, value) {

    if (ms(object))
        stop ("clustertrap requires single-session 'traps' object")

    if (!(inherits(object, 'traps')))
        stop ("requires clustered 'traps' object")

    if (length(value) > 0) {
        if (length(value) != nrow(object))
            stop ("requires one clustertrap per detector or detector vertex")
    }
    else
        value <- NULL

    structure (object, clustertrap = value)
}

'signalframe<-' <- function (object, value) {
    if (is.null(value)) {
        attr(object, 'signalframe') <- NULL
        object
    }
    else {
        value <- as.data.frame(value)
        if (nrow(value) != sum(abs(object)))
            stop ("requires one row per detection")
        if (!('signal' %in% names(value)))
            stop ("value does not contain column 'signal'")
        if (!(detector(traps(object)) %in% c('signal','signalnoise')) |
            !(inherits(object,'capthist')))
            stop ("requires 'capthist' object with 'signal' or 'signalnoise' detector")
        if (ms(object))
            stop ("requires single-session 'capthist' object")
        structure (object, signalframe = value)
    }
}

'signal<-' <- function (object, value) {
    if (length(value) != sum(abs(object)))
        stop ("requires one signal per detection")
    if (!(detector(traps(object)) %in% c('signal','signalnoise')) |
        !(inherits(object,'capthist')))
        stop ("requires 'capthist' object with 'signal' or 'signalnoise' detector")
    if (ms(object))
        stop ("requires single-session 'capthist' object")
    sf <- attr(object, 'signalframe')
    if (is.null(sf)) {
        sf <- data.frame(signal = value)
    }
    else
        sf$signal <- value
    structure (object, signalframe = sf)
}

'noise<-' <- function (object, value) {
    if (length(value) != sum(abs(object)))
        stop ("requires one noise per detection")
    if (!(detector(traps(object)) %in% c('signalnoise')) |
        !(inherits(object,'capthist')))
        stop ("requires 'capthist' object with 'signalnoise' detector")
    if (ms(object))
        stop ("requires single-session 'capthist' object")
    sf <- attr(object, 'signalframe')
    if (is.null(sf))
        sf <- data.frame(noise = value)
    else
        sf$noise <- value
    structure (object, signalframe = sf)
}

'traps<-' <- function (object, value) {
    if (!is(value,'traps'))
        stop ("'traps' object required for replacement")

    ## MODIFIED 2010 04 27
    if (ms(object)) {
        nsess <- length(object)
        temp <- vector(mode='list', nsess)
        if (nsess != length(value))
            stop ("replacement value has wrong length")
        for (i in 1:nsess) temp[[i]] <- `traps<-`(object[[i]], value[[i]])
        class(temp) <- c('list', 'traps')
        temp
    }
    else {
        structure (object, traps = value)
    }
}

'session<-' <- function (object, value) {
    if (ms(object)) {
       if (length(value) != length(object))
           stop ("invalid replacement value")
       for (i in 1:length(object)) session(object[[i]]) <- value[i]   ## 2010 03 26
       structure (object, names = as.character(value))
    }
    else {
        if (length(value) > 1)
            stop ("requires only one session name per session")
        structure (object, session = as.character(value))
    }
}
###############################################################################

######################################
## Class : traps
## defines an array of detectors
## detector may be 'single', 'multi', 'proximity' etc.
######################################

## 2011-10-10 make.grid, makepoly etc. moved to make.grid.R

rotate.traps <- function (object, degrees, centrexy=NULL, ...)
{
##    if (ms(object)) lapply(object, rotate.traps, degrees, centrexy, ...)
##    else

  rotatefn <- function (xy) {
    # about centre
    x <- xy[1] - centrexy[1]
    y <- xy[2] - centrexy[2]
    x2 <- x * cos(theta) + y * sin(theta) + centrexy[1]
    y2 <- - x * sin(theta) + y * cos(theta) + centrexy[2]
    c(x2,y2)
  }

  if (abs(degrees)>0) {
    if (is.null(centrexy)) centrexy <- c(mean(object$x), mean(object$y))
    theta <- 2*pi*degrees/360 # convert to radians

    traps2 <- data.frame(t(apply (object,1,rotatefn)))
    names(traps2) <- c('x','y')
    attr(traps2,'class')  <- c('traps', 'data.frame')
    detector(traps2)      <- detector(object)
    if (!is.null(usage(object)))
        usage(traps2)         <- usage(object)
    if (!is.null(covariates(object)))
        covariates(traps2) <- covariates(object)
    if (!is.null(timevaryingcov(object)))
        timevaryingcov(traps2) <- timevaryingcov(object)
    if (!is.null(polyID(object)))   ## includes transectID
        polyID(traps2)    <- polyID(object)
  }
  else traps2 <- object
  traps2
}
###############################################################################

shift.traps <- function (object, shiftxy, ...)
{
  object$x <- object$x + shiftxy[1]
  object$y <- object$y + shiftxy[2]
  object
}
###############################################################################

shift.mask <- function (object, shiftxy, ...)
{
  object$x <- object$x + shiftxy[1]
  object$y <- object$y + shiftxy[2]
  object
}
###############################################################################

flip.traps <- function (object, lr=F, tb=F, ...) {

##    if (ms(object)) lapply(object, flip.traps, lr, tb, ...)
##    else

    if (is.logical(lr)) {
        if (lr) object$x <- 2 * mean(object$x) - object$x  ## flip about centre
    } else
        if (is.numeric(lr)) object$x <- 2*lr - object$x  ## flip about lr

    if (is.logical(tb)) {
        if (tb) object$y <- 2 * mean(object$y) - object$y  ## flip about centre
    } else
        if (is.numeric(tb)) object$y <- 2*tb - object$y  ## flip about tb

    object
}
###############################################################################

rotate.popn <- function (object, degrees, centrexy=NULL, ...)
{
  rotatefn <- function (xy) {
    # about centre
    x <- xy[1] - centrexy[1]
    y <- xy[2] - centrexy[2]
    x2 <- x * cos(theta) + y * sin(theta) + centrexy[1]
    y2 <- - x * sin(theta) + y * cos(theta) + centrexy[2]
    c(x2,y2)
  }
  bbox <-  attr(object, 'boundingbox')
  if (abs(degrees)>0) {
    if (is.null(centrexy)) centrexy <- c(mean(bbox$x), mean(bbox$y))
    theta <- 2*pi*degrees/360 # convert to radians

    popn2 <- data.frame(t(apply (object,1,rotatefn)))

    object[,] <- popn2[,]
    attr(object, 'boundingbox') <- data.frame(rotate (bbox, degrees, centrexy))
  }
  object
}
###############################################################################

shift.popn <- function (object, shiftxy, ...)
{
  object$x <- object$x + shiftxy[1]
  object$y <- object$y + shiftxy[2]
  bbox <-  attr(object, 'boundingbox')
  attr(object, 'boundingbox') <- data.frame(shift(bbox, shiftxy))
  object
}
###############################################################################

flip.popn <- function (object, lr=F, tb=F, ...) {
    bbox <- attr(object, 'boundingbox')
    if (is.logical(lr)) {
        if (lr) object$x <- 2 * mean(bbox$x) - object$x  ## flip about centre
    } else
        if (is.numeric(lr)) object$x <- 2*lr - object$x  ## flip about lr

    if (is.logical(tb)) {
        if (tb) object$y <- 2 * mean(bbox$y) - object$y  ## flip about centre
    } else
        if (is.numeric(tb)) object$y <- 2*tb - object$y  ## flip about tb
    attr(object, 'boundingbox') <- data.frame(flip(bbox, lr=lr, tb=tb, ...))
    object
}
###############################################################################

print.traps <- function(x, ...) {
    if (ms(x)) {
        for (i in 1:length(x)) {
            cat('\n')
            if (!is.null(names(x)))
                cat(names(x)[i], '\n')
            print (x[[i]], ...)
        }
        invisible()
    }
    else {
        temp <- data.frame(row.names=attr(x,'row.names'), x=x$x, y=x$y)
        print(temp, ...)
    }
}
###############################################################################

subset.traps <- function (x, subset = NULL, occasions = NULL, ...) {
    # subset may be numeric index or logical

    if (ms(x)) {
        temp <- lapply(x, subset.traps, subset=subset, occasions = occasions, ...)
        class(temp) <- c('list', 'traps')
        temp
    }
    else {
        ## 2011-03-29
        if (is.null(subset)) {
            if (detector(x) %in% .localstuff$polydetectors)
                subset <- 1:length(levels(polyID(x)))
            else
                subset <- 1:nrow(x)
        }

        ## polygon & transect objects subset by whole polygons or transects
        ## 2011-01-24
        rowsubset <- subset  ## default
        if (detector(x) %in% c('polygon', 'polygonX', 'telemetry'))
            rowsubset <- as.character(polyID(x)) %in% levels(polyID(x))[subset]
        if (detector(x) %in% c('transect', 'transectX'))
            rowsubset <- as.character(transectID(x)) %in% levels(transectID(x))[subset]
        ## apply subsetting
        temp <- x[rowsubset,,drop=F]
        class(temp) <- c('traps','data.frame')
        detector(temp) <- detector(x)

        ## 2011-05-09
        if (!is.null(clusterID(x))) {
            clusterID(temp) <- factor(clusterID(x)[rowsubset])
            clustertrap(temp) <- factor(clustertrap(x)[rowsubset])
        }

        ## restore polyiD, transectID, usage, covariates
        if (detector(x) %in% c('polygon', 'polygonX', 'telemetry'))
            polyID(temp) <- factor(polyID(x)[rowsubset])
        if (detector(x) %in% c('transect', 'transectX'))
            transectID(temp) <- factor(transectID(x)[rowsubset])

        if (!is.null(usage(x))) {
            if (is.null(occasions))
                occasions <- 1:ncol(usage(x))
            usage(temp) <- usage(x)[subset,occasions,drop=F]
        }
        if (!is.null(markocc(x))) {
            if (is.null(occasions))
                occasions <- 1:length(markocc(x))
            markocc(temp) <- markocc(x)[occasions]
        }
        if (!is.null(covariates(x)))
            covariates(temp) <- covariates(x)[subset,,drop=F]

        if (!is.null(timevaryingcov(x))) {
            timevaryingcov(temp) <- lapply(timevaryingcov(x),
                function(y) y[occasions,drop=FALSE])
        }

        ## 2012-08-07
        if (detector(x) %in% .localstuff$pointdetectors) {
            spacing(temp) <- spacing(temp, recalculate = TRUE)
        }

        temp
    }
}
###############################################################################

rbind.traps <- function (..., renumber = TRUE, addusage ) {
# combine 2 or more traps objects
# what if multi-session?
    allargs <- list(...)
    check <- function (x) {
        result <- FALSE
        if (!is(x,'traps'))
            stop ("all arguments must be 'traps' objects")
        if (is.null(covariates(x)) != is.null(covariates(allargs[[1]]) ))
            stop ("covariates must be provided for all or none")
        if (detector(x) != detector(allargs[[1]]))
            warning ("detector types vary; using first")
        if (is.null(timevaryingcov(x)) != is.null(timevaryingcov(allargs[[1]]) ))
            result <- TRUE
        result
    }

    if (length(allargs) <= 1) {
        if (length(allargs)==1 & ms(allargs[[1]])) {
            allargs <- allargs[[1]]
            temp <- do.call(rbind.data.frame, allargs)
        }
        else
            stop ("requires more than one traps object")
    }
    else {
        tvc <- sapply (allargs, check)
        if (any(tvc))
            warning("time-varying covariates differ; using first")
        temp <- rbind.data.frame(...)
    }
    class(temp) <- c('traps', 'data.frame')

    ## 2010 07 05
    tempdet <- lapply(allargs, detector)
    if (length(unique(sapply(tempdet, detector))) >1 )
        ## stop ("cannot combine detector types - change with reduce()")
        ## 2011-09-26
        warning ("combining detector types - maybe change first with reduce()?")
    detector(temp) <- tempdet[[1]]

    ## 2011-02-07
    ## code for originating traps object
    ## uses 'polyID = transectID' behaviour of polyID()
    if (detector(allargs[[1]]) %in% .localstuff$polydetectors) {
        oldtrapsID <- rep(1:length(allargs), sapply(allargs, nrow))
        newpolyID <- unlist(lapply(allargs, function(x) as.character(polyID(x))))
        newpolyID <- factor(paste (oldtrapsID,newpolyID, sep='.'))
        polyID(temp) <- newpolyID
    }

    ## covariates 2010 07 05, 2010-08-28
    tempcov <- lapply(allargs, covariates)
    common <- Reduce(intersect, lapply(tempcov, names))
    tempcov <- lapply(tempcov, function(x) x[,common, drop = FALSE])
    covariates(temp) <- do.call(rbind, tempcov)
    timevaryingcov(temp) <- timevaryingcov(allargs[[1]])

    ## clusters  2011-04-12
    tempclus <- lapply(allargs, clusterID)
    if (!is.null(tempclus[[1]])) {
        tempclus <- lapply(tempclus, function(x) as.numeric(as.character(x)) )
        clusterID(temp) <- do.call(c, tempclus)
        clustertrap(temp) <- unlist(lapply(allargs, clustertrap))
    }
    else {
        clusterID(temp) <- NULL
        clustertrap(temp) <- NULL
    }

    ## usage
    ## tweaked 2013-01-16
    tempusage <- lapply(allargs, usage)
    if (any(!sapply(tempusage, is.null))) {
        nocc <- unique(unlist(sapply(tempusage, ncol)))
        nr <- unlist(sapply(allargs, nrow))
        if (length(nocc)>1)
            warning ("varying number of occasions; using maximum")
        nocc <- max(nocc)
        fillmissing <- function(x, nr) {
            flush.console()
            if(is.null(x))
                x <- matrix(1, nrow = nr, ncol = nocc)
            else {
                if (ncol(x) < nocc) {
                    x <- cbind(x, matrix(0, nrow = nr, ncol = nocc-ncol(x)))
                }
            }
            ## 2016-01-07 temporary fix
            ## dimnames(x)[[2]] <- 1:ncol(x)
            dimnames(x) <- list(NULL, 1:ncol(x))
            x
        }
        for (i in 1: length(tempusage)) {
            tempusage[[i]] <- fillmissing(tempusage[[i]], nr[i])
        }
        usage(temp) <- do.call(rbind, tempusage)
    }
    ## 2014-03-19
    else if (!missing(addusage)) {
        if (length(addusage) == 1)
            addusage <- rep(addusage, length(allargs))
        if ((length(addusage) != length(allargs)) | any(addusage<1))
            stop ("invalid occasion numbers in addusage")
        makeusage <- function(nr,nc) {
            tmp <- matrix(0, nrow = nr, ncol = nocc)
            tmp[,1:nc] <- 1
            tmp
        }
        nocc <- max(addusage)
        nt <- sapply(allargs, nrow)
        tempusage <- mapply(makeusage, nt, addusage)
        usage(temp) <- do.call(rbind, tempusage)
    }

    tempmarkocc <- lapply(allargs, markocc)
    mo <- !sapply(tempmarkocc, is.null)
    if (any(mo)) {
        if (!all(mo))
            warning ("using first non-null markocc")
        wmo <- which(mo)
        mo1 <- tempmarkocc[[wmo[1]]]
        if (length(wmo) > 1) {
            for (i in 2:length(wmo))
                if (any(tempmarkocc[[i]] != mo1))
                    stop ("conflicting markocc")
        }
        markocc(temp) <- mo1
    }

    tn <- sapply(allargs, row.names, simplify=F)


 #   if (length(unique(unlist(tn))) != length(unlist(tn))) {  # renumber
        if (renumber) {
            if (detector(temp) %in% .localstuff$polydetectors) {
                if (!is.null(polyID(temp)))  ## polyID also stores transectID
                    polyID(temp) <- factor(as.numeric(polyID(temp)))
                temp <- renamepolyrows(temp)    ## see read.traps
            }
            else
                row.names(temp) <- 1:nrow(temp)
        }
        else {
            for (i in 1:length(tn)) tn[[i]] <- paste(tn[[i]],i,sep='.')
            row.names(temp) <- unlist(tn)
        }
#    }
    if (!is.null(usage(temp)))
        if (nrow(usage(temp))>0)
            row.names(usage(temp)) <- row.names(temp)

    if (!is.null(covariates(temp))) {
        if (nrow(covariates(temp))>0)
            row.names(covariates(temp)) <- row.names(temp)
    }

    ## 2012-08-07, 2012-08-24
    if (detector(temp) %in% .localstuff$pointdetectors) {
        spacing(temp) <- spacing(temp, recalculate = TRUE)
    }

    temp
}
###############################################################################

plot.traps <- function(x,
    border = 100,
    label = FALSE,
    offset = c(6,6),
    add = FALSE,
    hidetr = FALSE,
    detpar=list(),
    txtpar=list(),
    bg = 'white',
    gridlines = TRUE,
    gridspace = 100,
    gridcol = 'grey',
    markused = FALSE,
    markvarying = FALSE,
    markvertices = FALSE,
    labelclusters = FALSE,
    ... )
{
#### NEED TO HANDLE CLUSTER, CLUSTERTRAP 2011-04-12

    if (ms(x)) {
        lapply(x, plot.traps,
            border, label, offset,
            add, hidetr, detpar, txtpar,
            bg, gridlines, gridspace, gridcol,
            markused, markvarying, markvertices, labelclusters, ...)
    }
    else {
        buff <- c(-border,+border)
        offsety <- ifelse (length(offset)==2, offset[2], offset[1])
        dcol <- 'red'
        detpar <- replacedefaults (list(col=dcol, pch=3, cex=0.8), detpar)
        txtpar <- replacedefaults (list(col='blue', cex=0.7), txtpar)
        if (is.logical(markvertices))
            markvertices <- markvertices * 2  ## 0 or 2

        if (!is.null(usage(x))) {
            used <- apply(attr(x,'usage'),1,function(z) any(z>0))
            varying <- used * apply(attr(x,'usage'),1,function(z) any(z==0))
        }
        else {
            used  <- rep(TRUE, nrow(x))
            varying <- rep(FALSE, nrow(x))
        }
        initialpar <- par(detpar)
        on.exit(par(initialpar))

        if (!add) {
            par(bg=bg)
            ## axes = FALSE blocks bty = 'o' 2011-05-08
            eqscplot (x$x, x$y, xlim=range(x$x)+buff, ylim=range(x$y)+buff,
                xlab='', ylab='', type='n', axes=F, ...)
            if (gridlines) {
                xl <- range(x$x)+buff
                yl <- range(x$y)+buff
                strtx <- gridspace * floor(xl[1]/gridspace)
                strty <- gridspace * floor(yl[1]/gridspace)
                finx  <- gridspace * (floor(xl[2]/gridspace) + 1)
                finy  <- gridspace * (floor(yl[2]/gridspace) + 1)
                for (xi in seq(strtx, finx, gridspace))
                    segments(xi, strty, xi, finy, col=gridcol)
                for (yi in seq(strty, finy, gridspace))
                    segments(strtx, yi, finx, yi, col=gridcol)
            }
        }
        plotvertices <- function (df) {
            if (markvertices == 1)
                i <- c(1, nrow(df))
            else
                i <- 1:nrow(df)
            points(df$x[i], df$y[i], pch = detpar$pch, bg = bg, col=detpar$col)
        }

        if (hidetr==F) {
            if (detector(x) %in% c('polygon','polygonX', 'telemetry')) {
                templist <- split (x, levels(polyID(x)), prefix='')
                lapply(templist, function (y)
                       polygon (y$x, y$y, col=detpar$col, density=0))
                if (markvertices > 0) {
                    lapply(templist, plotvertices)
                }
                if (label) for (k in 1:length(templist)) {
                    if (detector(x) %in% c('polygon','polygonX')) {
                        msk <- suppressWarnings(make.mask(templist[[k]], buffer = 0, poly = templist[[k]], nx = 32))
                        xbar <- mean(msk$x)
                        ybar <- mean(msk$y)
                    }
                    else {
                        xbar <- mean(range(templist[[k]]$x))
                        ybar <- mean(range(templist[[k]]$y))
                    }
                    text (xbar+offset[1], ybar+offsety, names(templist)[k])
                }
            }
            else if (detector(x) %in% c('transect','transectX')) {
                templist <- split (x, levels(transectID(x)), prefix='')
                lapply(templist, function (df) lines (df$x, df$y, col=detpar$col))
                if (markvertices > 0) {
                    lapply(templist, plotvertices)
                }
                if (label) for (k in 1:length(templist)) {
                    xbar <- mean(range(templist[[k]]$x))
                    ybar <- mean(range(templist[[k]]$y))
                    text (xbar+offset[1], ybar+offsety, names(templist)[k])
                }
            }
            else {
                points (x$x, x$y)
                if (markused) {
                    points (x$x[used], x$y[used], pch = 1, cex = 0.8)
                }
                if (markvarying & any(varying)) {
                    points (x$x[varying], x$y[varying], pch = 16, cex = 0.8)
                }
            }
            par(txtpar)
            if (label && !(detector(x) %in% .localstuff$polydetectors)) {
                text (x$x+offset[1], x$y+offsety, rownames(x))
            }
            if (labelclusters && !(detector(x) %in% .localstuff$polydetectors)) {
                if (is.null(clusterID(x)) | is.null(clustertrap(x)))
                    stop ("require clustered traps to label with clusterID")
                cl1 <- clustertrap(x) == 1
                text (x$x[cl1]+offset[1], x$y[cl1]+offsety, clusterID(x)[cl1])
            }
            #par(initialpar)   # restore
        }
        invisible()
    }
}
###############################################################################

summary.traps <- function(object, getspacing = TRUE, ...) {

#### NEED TO HANDLE CLUSTER, CLUSTERTRAP 2011-04-12
    if (ms(object)) lapply(object, summary.traps, getspacing = getspacing, ...)
    else {
        if (is.null(object$x) | is.null(object$y))
            stop ("not a valid 'traps' object")
        nd <- nrow(object)
        np <- NA
        sumusage <- NULL
        rangeusage <- NULL
        sumcovar <- NULL

        if (detector(object) %in% .localstuff$polydetectors) {
            spacing <- NA
        }
        else {
            spacing <- spacing(object, recalculate = getspacing)

            # unclear why following is not for polydetectors?  FAULTY 2012-12-18
            if (is.factor(covariates(object))) {
                susage <- by(usage(object), covariates(object), function(y) apply(y,2,sum))
                sumusage <- matrix(unlist(susage), byrow = T, nrow = length(susage))
                dimnames(sumusage) <- list(levels(covariates(object)), names(susage[[1]]))
            }
            else if (!is.null(usage(object)))  {
                sumusage <- apply(usage(object)>0, 2, sum)
                nocc <- ncol(usage(object))
                rangeusage <- matrix(apply(usage(object), 2, range), ncol=nocc)
                dimnames(rangeusage) <- list(c('min','max'), 1:nocc)
            }
            tempcovar <- covariates(object)

            if (!is.null(tempcovar))
                if ((nrow(tempcovar)>0) & (ncol(tempcovar)>0)) ## amended to check ncol 2009 09 18
                    sumcovar <- summary(tempcovar, ...)
        }

        ## defaults
        area <- NA
        totallength <- NA
        np = ndetector(object)
        spacex <- NA
        spacey <- NA
        ## unblock these 2010-12-17
        sumusage <- NULL
        sumcovar <- NULL
        xrange <- range(object$x)
        yrange <- range(object$y)
        if (detector(object) %in% c('polygon', 'polygonX', 'telemetry'))
            area <- searcharea(object)
        if (detector(object) %in% c('transect', 'transectX'))
            totallength <- transectlength(object)
        tvc <- timevaryingcov(object)

        temp <- list (
          detector = detector(object),
          ndetector = nd,
          npart = np,
          xrange = xrange,
          yrange = yrange,
          spacing = spacing,
          area  = area,
          totallength = totallength,
          usage = sumusage,
          markocc = markocc(object),
          rangeusage = rangeusage,
          covar = sumcovar,
          tvc = tvc
        )

        class(temp) <- 'summary.traps'
        temp
    }
}
###############################################################################

print.summary.traps <- function (x, terse = FALSE, ...) {

#### NEED TO HANDLE CLUSTER, CLUSTERTRAP 2011-04-12

    if (!terse)
    cat ('Object class     ', 'traps', '\n')
    cat ('Detector type    ', x$detector, '\n')
    if (x$detector %in% c('polygon', 'polygonX', 'telemetry')) {
        cat ('Number vertices  ', x$ndetector-x$npart, '\n')  ## assume each polygon closed
        cat ('Number polygons  ', x$npart, '\n')
        cat ('Total area       ', sum(x$area), 'ha \n')
    }
    else if (x$detector %in% c('transect', 'transectX')) {
        cat ('Number vertices  ', x$ndetector, '\n')
        cat ('Number transects ', x$npart, '\n')
        cat ('Total length     ', x$totallength, 'm \n')
    }
    else {
        cat ('Detector number  ', x$ndetector, '\n')
        cat ('Average spacing  ', x$spacing, 'm \n')
    }
    cat ('x-range          ', x$xrange, 'm \n')
    cat ('y-range          ', x$yrange, 'm \n')
    if (!is.null(x$rangeusage)) {
        cat ('\n Usage range by occasion\n')
        print(x$rangeusage, ...)
    }
    if (!is.null(x$markocc)) {
        cat ('\nMarking occasions\n')
        mo <- matrix(as.numeric(x$markocc), byrow=T, nrow = 1, 
                     dimnames = list('',1:length(x$markocc)))
        print(mo, ...)
    }
    if (!terse) {
        if (!is.null(x$covar)) {
            cat ('\n')
            cat ('Summary of covariates', '\n')
            print(x$covar, ...)
        }
        if (!is.null(x$usage)) {
            cat ('Usage>0 by occasion', '\n')
            print(x$usage, ...)
        }
        if (!is.null(x$tvc)) {
            cat ('Time-varying covariate(s) name : columns', '\n')
            for (i in 1:length(x$tvc)) cat(names(x$tvc)[[i]], ':', x$tvc[[i]], '\n')
        }
    }
}

###############################################################################

####################################
## Class : capthist
## capture data
####################################

###############################################################################

plot.popn <- function (x, add = FALSE, frame = TRUE, circles = NULL, ...) {

    if (ms(x)) {
        ## force shared frame
        temp <- do.call(rbind, lapply(x, function(y) attr(y,'boundingbox')))
        vertices <- apply(temp,2,range)
        for (i in 1:length(x)) attr(x,'boundingbox') <- vertices
        lapply (x, plot, add, frame, circles, ...)
        invisible()
    }
    else {
        vertices <- attr(x,'boundingbox')

        if (add==FALSE)
        {
            if (frame)
                eqscplot (x$x, x$y, xlab='', ylab='', xlim=range(vertices$x),
                    ylim=range(vertices$y), type='n', axes = FALSE, ...)
            else
                eqscplot (x$x, x$y, xlab='', ylab='', type='n', axes = FALSE,
                          ...)
        }
        if (is.null(circles) | (nrow(x) == 0))    ## second condition 2011-09-14
            points (x$x, x$y, ...)
        else {
            if (length(circles) == 1)
                circles <- rep(circles, nrow(x))
            symbols (x$x, x$y, circles = circles, inches = FALSE,
                add = TRUE, ...)
        }
        if (frame) {
            if (!is.null(attr(x,'polygon')))
                polygon (attr(x,'polygon'))
            else
                polygon (vertices)
        }
    }
}

###############################################################################
## 2015-12-18 deparse.level added

rbind.popn <- function (..., renumber = TRUE) {
## combine 2 or more popn objects
## ... may be a single list object

    dots <- match.call(expand.dots = FALSE)$...
    allargs <- list(...)
    
    ## ad hoc fix for occasional inclusion of deparse.level 2015-12-18
    if (length(dots)>length(allargs))
        dots <- dots[-1]
    
    names(allargs) <- lapply(dots, as.character)

    if ((length(dots)==1) & (!inherits(allargs[[1]],'popn'))) allargs <- allargs[[1]]

    if (length(allargs)==1) return(allargs[[1]])

    ## check input
    check <- function (x) {
        if (!is(x,'popn'))
            stop ("all arguments must be 'popn' objects")
        if (is.null(covariates(x)) != is.null(covariates(allargs[[1]]) ))
            stop ("covariates must be provided for all or none")
    }
    sapply (allargs, check)

    ## row names
    an <- unlist(sapply(allargs, row.names, simplify=F))
    names(an) <- NULL
    if (any(duplicated(an))) # renumber
    {
        if (renumber) rn <- 1:length(an)
        else {
            for (i in 1:length(an)) an[[i]] <- paste(an[[i]],i,sep='.')
            rn <- unlist(an)
        }
    }
    else rn <- an

    ## construct output
    animals <- data.frame(abind (allargs, along=1), row.names = rn)
    names(animals) <- c('x','y')
    class(animals) <- c('popn', 'data.frame')
    ## following 2 lines modified 2010-06-13
    attr(animals, 'Ndist') <- 'user'
    attr(animals, 'model2D') <- attr(allargs[[1]], 'model2D')
    xl <- range(sapply(allargs, function(x) attr(x,'boundingbox')$x))
    yl <- range(sapply(allargs, function(x) attr(x,'boundingbox')$y))
    attr(animals, 'boundingbox') <- expand.grid (x=xl,y=yl)[c(1,3,4,2),]
    if (!is.null(covariates(allargs[[1]]))) {
        cov <- lapply(allargs, function(x) covariates(x))
        covariates(animals) <- data.frame(abind(cov, along=1), row.names=rn)
    }
    animals
}
###############################################################################

subset.popn <- function (x, subset = NULL, sessions = NULL, poly = NULL,
    poly.habitat = TRUE, keep.poly = TRUE, renumber = FALSE, ...)
## x - popn object
## subset - vector (character, logical or integer) to subscript rows (dim(1)) of x
## sessions - vector (integer or logical) to subscript sessions
{
    if (ms(x)) {
        if (!is.null(sessions))
            x <- x[sessions]
        out <- vector('list')
        for (i in 1:length(x)) {
            if (is.list(subset))
                sset <- subset[[i]]
            else
                sset <- subset
            if (is.null(subset))  ## 2015-03-17
                out[[i]] <- x[[i]]
            else
                out[[i]] <- subset(x[[i]], subset[[i]], NULL, renumber, ...)
        }
        class(out) <- c('list','popn')
        out
    }
    else {
        #-------------------------
        # default subset is all
        if (is.null(subset))
            subset <- 1:nrow(x)
        #-------------------------
        # restrict to a polygon
        # added 2011-10-20

        if (!is.null(poly)) {
            OK <- pointsInPolygon(x, poly)
            if (!poly.habitat)
                OK <- !OK
            subset <- subset[OK]
        }
        #-------------------------
        # apply subsetting
        pop <- x[subset,]
        #-------------------------
        if (renumber)
            rownames(pop) <- 1 : nrow(pop)

        class(pop) <- c('popn', 'data.frame')
        attr(pop, 'Ndist') <- NULL     ## no longer known
        attr(pop, 'model2D') <- NULL   ## no longer known
        attr(pop, 'boundingbox') <- attr(x, 'boundingbox')
        if (!is.null(poly) & keep.poly) {
            attr(pop, 'polygon') <- poly
            attr(pop, 'poly.habitat') <- poly.habitat
        }
        if (!is.null(covariates(x))) {
            covariates(pop) <- covariates(x)[subset,,drop=FALSE]
        }
        pop
    }
}

###############################################################################

subset.capthist <- function (x, subset=NULL, occasions=NULL, traps=NULL,
    sessions=NULL, cutval=NULL, dropnullCH=TRUE, dropnullocc=FALSE,
    dropunused = TRUE, droplowsignals = TRUE, dropNAsignals = FALSE,
    cutabssignal = TRUE, renumber=FALSE, ...)  {

## modified 2011-11-14; 2012-02-11

## x - capthist object (array with 2 or 3 dimensions)
## subset - vector (character, logical or integer) to subscript rows (dim(1)) of x
## occasions - vector (integer or logical) to subscript occasions
## traps - vector (character, integer or logical) to subscript rows of traps object
## sessions - vector (integer or logical) to subscript sessions
    if (ms(x)) {
        if (is.null(sessions)) sessions <- 1:length(x)
        temp <- lapply (x[sessions], subset,
            subset = subset,
            occasions = occasions,
            traps = traps,
            sessions = sessions,   ## inserted 2009 10 01
            cutval = cutval,
            dropnullCH = dropnullCH,
            dropnullocc = dropnullocc,
            dropunused = dropunused,
            droplowsignals = droplowsignals,
            renumber = renumber, ...)
        class(temp) <- c('list', 'capthist')
        if (length(temp) == 1) temp <- temp[[1]]  ## 2009 09 25
        return(temp)
    }
    else {
    
        detector <- detector(traps(x))
        dim3 <- length(dim(x)) == 3
        nk <- ndetector (traps(x))
        if (is.logical(subset) & (length(subset) != nrow(x)))
            stop ("if 'subset' is logical its length must match number of animals")
        if (is.logical(occasions) & (length(occasions) != ncol(x)))
            stop ("if 'occasions' is logical its length must match number of occasions")
        if (is.logical(traps) & (length(traps) != nk))
            stop ("if 'traps' is logical its length must match number of detectors")
        if (is.null(occasions)) occasions <- 1:ncol(x)
        if (is.null(traps))  traps <- 1:nk
        if (is.null(subset)) subset <- 1:nrow(x)

        #############################################
        ## coerce subset, traps, occasions to logical
        if (is.character(subset))
            subset <- dimnames(x)[[1]] %in% subset
        else
            if (!is.logical(subset))
                subset <- (1:nrow(x)) %in% subset
        if (is.character(traps))
            traps <- rownames(traps(x)) %in% traps
        else
            if (!is.logical(traps))
                traps <- (1:nk) %in% traps
        if (!is.logical(occasions))
            occasions <- (1:ncol(x)) %in% occasions

        #####################################
        ## signaldf is used later...
        if (detector %in% c('signal','signalnoise')) {
            signaldf <- data.frame(trap = trap(x, names = F),
                                   occ = occasion(x),
                                   ID = animalID(x, names = F),
                                   attr(x, 'signalframe'))
        }

        #####################################
        ## apply signal threshold if relevant
        if (detector %in% c('signal','signalnoise')) {
            if (is.null(cutval)) cutval <- attr(x, 'cutval')
            if (cutval < attr(x, 'cutval'))
                stop ("cannot decrease 'cutval'")
            if (cutabssignal) {
                signalOK <- (signal(x) >= cutval)
            }
            else {
                if (is.null(noise(x)))
                    stop("could not find noise for relative signal cut")
                signalOK <- ((signal(x)-noise(x)) >= cutval)
            }
            signalOK <- ifelse(is.na(signalOK), !dropNAsignals, signalOK)
            newcount <- table(
                factor(animalID(x, names = FALSE), levels = 1:nrow(x))[signalOK],
                factor(occasion(x), levels = 1:ncol(x))[signalOK],
                factor(trap(x, names = FALSE), levels = 1:nk)[signalOK])
            if (nrow(x)>0)  ## 2013-08-17
            x[] <- newcount * sign(x)  ## retain deads, in principle
        }

        ###########################
        ## condition missing values
        x[is.na(x)] <- 0
        #################################
        ## optionally drop traps never used on the specified occasions
        if (dropunused && !is.null(usage(traps(x)))) {
            used <- apply(usage(traps(x))[,occasions, drop=F],1,sum) > 0
            traps <- traps & used
        }

        ##############################################
        ## which cells are retained?
        ## for signalframe 2012-02-11
        if (dim3) {
            i <- x
            if (nrow(x) > 0) {
                i[] <- 1:length(i)
                i <- i[subset, occasions, traps, drop = FALSE]
            }
        }

        #################################
        ## prepare to drop null histories
        if (dropnullCH) {
            if (dim3) {
                nonnull <- apply(abs(x[,occasions, traps, drop=FALSE]),1,sum) > 0
            }
            else {
                x[!(abs(x) %in% (1:nk)[traps])] <- 0
                nonnull <- apply(abs(x[,occasions, drop=F]),1,sum) > 0
            }
        }
        else
            nonnull <- rep(TRUE, nrow(x))
        subset <- subset & !is.na(subset) & nonnull
        if (nrow(x)==0) subset <- 0

        #################################
        ## perform main subset operation
        if (dim3) {
            temp <- x[subset, occasions, traps, drop = FALSE]
            ## parallel operation to track indices for dim3 case
        }
        else {
            temp <- x[subset, occasions, drop = FALSE]
            ## drop rejected trap sites; 'abs' added 2011-11-14
            temp[!(abs(temp) %in% (1:nk)[traps])] <- 0
        }

        nocc <- ncol(temp)

        #################################
        ## drop null occasions
        OK2 <- rep(T,nocc)
        if (nrow(temp)>0)
        if ((!(detector %in% c('signal','signalnoise'))) &&
            any( apply(abs(temp),2,sum) ==0)) {
            if (dropnullocc) {
                OK2 <- apply(abs(temp),2,sum) > 0
                if (dim3) {
                    temp <- temp[,OK2,, drop = FALSE]
                    i <- i[,OK2,, drop = FALSE]
                }
                else
                    temp <- temp[,OK2, drop = FALSE]
            }
            else
                warning ("no detections on occasion(s) ",
                    paste((1:nocc)[apply(abs(temp),2,sum) ==0], collapse=', '), "\n")
            nocc <- dim(temp)[2]  # refresh
        }
        ###################################
        ## attributes
        class(temp) <- 'capthist'
        traps(temp) <- subset (traps(x), traps)
        covariates(temp) <- covariates(x)[subset,,drop = FALSE]
        usage(traps(temp)) <- usage(traps(x))[traps, occasions,
            drop = FALSE][,OK2, drop = FALSE]  ## drop null occasions
        session(temp) <- session(x)
        attr(temp, 'n.mash') <- attr(x, 'n.mash')
        attr(temp, 'centres') <- attr(x, 'centres')

        ###################################
        ## mark-resight 2015-10-10,12
        if (!is.null(markocc(traps(x)))) {
            markocc(traps(temp)) <- markocc(traps(x))[occasions]
            if (!is.null(Tu(x))) {
                if (is.matrix(Tu(x))) {
                    Tu(temp) <- Tu(x)[traps, occasions, drop = FALSE]
                }
                else {
                    Tu(temp) <- Tu(x)
                    warning ('sightings not separated by occasion; all included')
                }
            }
            if (!is.null(Tm(x))) {
                if (is.matrix(Tm(x))) {
                    Tm(temp) <- Tm(x)[traps, occasions, drop = FALSE]
                }
                else {
                    Tm(temp) <- Tm(x)
                    warning ('nonID sightings not separated by occasion; all included')
                }
            }
        }
        
        ###################################
        ## 2012-10-20 for telemetry
        xylist <- telemetryxy(x)   ## attr(x, 'xylist')
        if (!is.null(xylist)) {
            xylist <- xylist[names(xylist) %in% rownames(x)[subset]]
            attr(temp, 'xylist') <- xylist
        }

        ###################################
        ## subset signal of signal capthist
        ## revised 2012-07-28
        if (detector %in% c('signal','signalnoise')) {
            if (nrow(x)>0) {
                if (!droplowsignals) {
                    if (any(x<=0))
                        stop ("droplowsignals = FALSE cannot be applied to CH",
                              " objects with incomplete detection")
                    signalOK <- TRUE
                }
            ## otherwise signalOK remains a logical vector with length equal
            ## to the original number of detections
            ## subset, traps and occasions are logical vectors 2011-01-21, 2011-11-14
            OK <- occasions[signaldf$occ] & subset[signaldf$ID] & traps[signaldf$trap] & signalOK
            signaldf <- signaldf[OK,, drop = FALSE]
            signaldf <- signaldf[order(signaldf$trap, signaldf$occ, signaldf$ID),, drop = FALSE]
            attr(temp, 'signalframe') <- signaldf[,-(1:3), drop = FALSE]
        }
            attr(temp, 'cutval') <- cutval
        }

        ############################################
        ## subset xy of polygon or transect capthist
        if (!is.null(xy(x))) {
            df <- data.frame(trap=trap(x, names = F),
                             occ=occasion(x),
                             ID=animalID(x,names = F),
                             x=xy(x)[,1], y=xy(x)[,2])
            ## subset, traps and occasions are logical vectors 2011-01-21, 2011-11-14
            OK <- occasions[df$occ] & subset[df$ID] & traps[df$trap]
            df <- df[OK,, drop = FALSE]
            df <- df[order(df$trap, df$occ, df$ID),, drop=FALSE]
            attr(temp, 'detectedXY') <- df[,c('x','y')]
        }

        ################################################
        # 2011-03-29, 2011-11-14
        # reassign trap numbers to allow for dropunused
        if (!dim3) {
            oldtrapnum <- 1:nk
            newtrapnum <- match(oldtrapnum, oldtrapnum[traps])
            trapsites <- as.numeric(abs(temp))
            temp[temp!=0] <- as.numeric(sign(temp[temp!=0])) * newtrapnum[trapsites]
        }

        ## renumber if desired
        if (renumber) {
            if (length(dim(temp))==3)
                dimnames(temp) <- list(1:nrow(temp),1:nocc,NULL)   # renew numbering
            else
                dimnames(temp) <- list(1:nrow(temp),1:nocc)
        }
        temp
    }
}

###############################################################################
##
## MS.capthist and rbind.capthist removed to rbind.capthist.R 13/9/2011
##
###############################################################################

sort.capthist <- function (x, decreasing = FALSE, by = '', byrowname = TRUE, ...) {
    if (ms(x)) {
        newx <- vector(mode='list')
        for (i in 1:length(x)) {
            xi <- subset(x, session=i)
            newx[[i]] <- sort(xi, decreasing=decreasing, by=by)
        }
        temp <- do.call(MS.capthist, newx)
        names(temp) <- names(x)
        session(temp) <- session(x)
        temp
    }
    else {
        if (is.character(by)) {
            if (by == '') {
                by <- vector(mode='list')   ## zero-length
            }
            else {
                if (!all(by %in% names(covariates(x))))
                    stop ("unrecognised sort field(s)")
                by <- as.list(covariates(x)[,by, drop=FALSE])
            }
        }
        else {
            by <- as.list(data.frame(by))
        }
        if (byrowname) by$rownames <- row.names(x)
        by$decreasing <- decreasing
        roworder <- do.call(order, by)
        newrownames <- rownames(x)[roworder]
        tempx <- x   ## for detectionorder
        tempx[tempx!=0] <- 1:length(animalID(x))
        if (length(dim(x))==2) {
            tempx[,] <- tempx[roworder,]
            x[,] <- x[roworder,]
        }
        else {
            tempx[,,] <- tempx[roworder,,]
            x[,,] <- x[roworder,,]
        }
        rownames(x) <- newrownames
        ## covariates
        if (!is.null(covariates(x))) {
            temp <-covariates(x)[roworder,,drop=F]
            rownames(temp) <- newrownames
            covariates(x) <- temp
        }
        
        detectionorder <- tempx[tempx!=0]
        ## signal
        sf <- attr(x, 'signalframe')
        if (!is.null(sf))
            attr(x, 'signalframe') <- sf[detectionorder,,drop=FALSE]
        ## xy
        if (!is.null(xy(x)))
            xy(x) <- xy(x)[detectionorder,,drop=FALSE]
        ## return
        x
    }
}

sort.mask <- function (x, decreasing = FALSE, by = '', byrowname = TRUE, ...) {
    if (ms(x)) {
        oldclass <- class(x)
        if (inherits(by, 'list'))
            temp <- mapply(sort, x, by, MoreArgs = list(decreasing = decreasing, ...), simplify = FALSE)  ## need to evaluate ...!!!
        else
            temp <- lapply(x, sort, decreasing = decreasing, by = by, ...)
        class (temp) <- oldclass
        temp
    }
    else {
        if (is.character(by)) {
            if (by == '') {
                by <- vector(mode='list')   ## zero-length
            }
            else {
                if (!all(by %in% names(covariates(x))))
                    stop ("unrecognised sort field(s)")
                by <- as.list(covariates(x)[,by, drop=FALSE])
            }
        }
        else {
            by <- as.list(data.frame(by))
        }
        if (byrowname) by$rownames <- row.names(x)
        by$decreasing <- decreasing
        roworder <- do.call(order, by)
        newrownames <- rownames(x)[roworder]
        temp <- x[roworder,]
        rownames(temp) <- newrownames
        ## covariates
        if (!is.null(covariates(x))) {
            tempcov <-covariates(x)[roworder,,drop=F]
            rownames(tempcov) <- newrownames
            covariates(temp) <- tempcov
        }
        
        attr(temp,'type')    <- 'sort'
        attr(temp,'meanSD')  <- attr(x,'meanSD')
        attr(temp,'area')    <- attr(x, 'area')
        attr(temp,'SLDF')    <- attr(x, 'SLDF')
        attr(temp, 'graph')  <- attr(x, 'graph')
        attr(temp,'spacing') <- attr(x, 'spacing')
        attr(temp,'boundingbox') <- attr(x,'boundingbox')
        class(temp) <- class(x)
        temp
    }
}

print.capthist <- function (x,..., condense = FALSE, sortrows = FALSE)
{
    ## recursive if list of capthist
    if (ms(x)) lapply (x, print.capthist, ..., condense = condense,
        sortrows = sortrows)
    else { # strip attributes, but why bother?
        cat('Session = ', session(x), '\n')
        if (is.null(detector(traps(x)))) {
            print.default(x, ...)
            invisible (x)
        }
        else
        if (detector(traps(x)) == 'unmarked') {
            temp <- apply(x>0, 2:3, sum)
            colnames(temp) <- rownames(traps(x))
            print.default(temp, ...)
            invisible (temp)
        }
        else
        if (detector(traps(x)) == 'presence') {
            ## temp <- apply(apply(x>0, 2:3, sum, drop = FALSE) > 0,2,any, drop = FALSE)*1.0
            ## 2015-05-10 retain occasion-specific data
            temp <- (apply(x>0, 2:3, sum, drop = FALSE) > 0)*1.0
            print.default(temp, ...)
            invisible (temp)
        }
        else
        if (condense & (detector(traps(x)) %in% c('proximity', 'count',
                                                  'signal','signalnoise'))) {
            temp <- apply(x, 3, function(y) y[apply(abs(y),1,sum)>0,, drop=F])
            trapnames <- rownames(traps(x))
            traps <- trapnames[rep(1:length(temp), sapply(temp, nrow))]
            detections <- as.matrix(abind(temp, along=1))
            temp <- data.frame(AnimalID = rownames(detections), Trap = traps, detections,
                stringsAsFactors = FALSE, row.names=NULL)
            names(temp)[3:ncol(temp)] <- 1:(ncol(temp)-2)
            if (sortrows) {
                lab <- temp$AnimalID
                if (suppressWarnings( all(!is.na(as.numeric(lab)))))
                    lab <- as.numeric(lab)
                temp <- temp[order(lab, temp$Trap ), ]
            }
            rownames(temp) <- 1:nrow(temp)
            print(temp, ...)
            invisible (temp)
        }
        else {
            temp <- array (x, dim=dim(x))
            dimnames(temp) <- dimnames(x)
            if (sortrows) {
                lab <- dimnames(temp)[[1]]
                if (suppressWarnings( all(!is.na(as.numeric(lab)))))
                    lab <- as.numeric(lab)
                if (length(dim(x)) == 3) temp <- temp[order(lab),,]
                else temp <- temp[order(lab),]
            }
            print.default(temp,...)
            invisible(temp)
        }
    }
}
############################################################################################

summary.capthist <- function(object, terse = FALSE, ...) {

    ## recursive if list of capthist
    if (ms(object)) {
        if (terse) {
            n     <- sapply(object, nrow)        # number caught
            nocc  <- sapply(object, ncol)        # number occasions
            ncapt <- sapply(object, function (xx)
                            if (detector(traps(xx)) %in% .localstuff$countdetectors)
                            sum(abs(xx))
                            else
                            sum(abs(xx)>0)
                            )
            ndet  <- sapply(traps(object), ndetector) # number traps
            temp  <- as.data.frame(rbind(nocc, ncapt, n, ndet))
            names(temp) <- names(object)
            rownames(temp) <- c('Occasions','Detections','Animals','Detectors')
            temp
        }
        else
            lapply (object, summary.capthist, ...)
    }
    else {

        traps <- traps(object)
        detector <- detector(traps)
        cutval <- attr(object, 'cutval')   # signal strength only

        nd <- length(traps$x)

        # ni, ui, fi, M, losses etc.
        nocc <- ncol(object)
        counts <- matrix(0, nrow = 8, ncol = nocc)
        signalsummary <- NULL
        if (nrow(object) > 0) {
            if (length(dim(object)) > 2) {
                tempx <- apply( object[,,,drop=F], c(1,2), function(x) sum(abs(x))>0)
                ## 2014-05-31
                tempx3 <-  apply( object[,,,drop=F], c(1,2), function(x) any(x<0))
                if (nocc>1) {  # distinction may not be needed...
                    counts [1,] <- apply(tempx, 2, function(x) sum(abs(x)>0) )
                    tempx2 <- apply(tempx, 1, function(x) cumsum(abs(x))>0)
                    counts [4,] <- apply(tempx2,1,sum)
                    counts [2,] <- c(counts[4,1],diff(counts[4,]))
                    counts [3,] <- tabulate(apply(tempx,1, function(x) sum(abs(x)>0)),nbins = nocc)
                    ## replaced counts [5,] <- apply(tempx,2, function(x) sum(x<0))
                    ## 2014-05-31
                    counts [5,] <- apply(tempx3, 2, sum)
                }
                else {
                    counts [1,1] <- sum(abs(tempx)>0)
                    counts [2,1] <- counts [1,1]
                    counts [3,1] <- counts [1,1]
                    counts [4,1] <- counts [1,1]
                    counts [5,1] <- sum(tempx<0)
                }
                if (detector %in% c('proximity','signal','signalnoise'))
                   counts [6,] <- apply(object[,,, drop=F], 2, function(x) sum(abs(x)>0))
                ## abs (x) added following 2 statements 2011-02-09
                if (detector %in% .localstuff$countdetectors) #c('count', 'polygon','transect','unmarked','presence'))
                   counts [6,] <- apply(object[,,, drop=F], 2, function(x) sum(abs(x)))
                tempt <- apply(object[,,,drop=F],c(2,3), function(x) sum(abs(x))>0)
                counts [7,] <- apply(tempt,1,sum)
            }
            else  {     ## .localstuff$exclusivedetectors
                if (nocc == 1) {
                    counts [1,1] <- sum(abs(object)>0)
                    counts [2,1] <- counts [1,1]
                    counts [3,1] <- counts [1,1]
                    counts [4,1] <- counts [1,1]
                    counts [5,1] <- sum(object<0)
                }
                else {
                    counts [1,] <- apply(object[,,drop=F], 2, function(x) sum(abs(x)>0))
                    if (nrow(object) > 0) {
                        tempM <- apply(object[,,drop=F], 1, function(x) cumsum(abs(x))>0)
                        counts [4,] <- apply(tempM,1,sum)
                    }
                    counts [2,] <- c(counts[4,1], diff(counts[4,]))
                    counts [3,] <- tabulate(apply(object[,,drop=F],1, function(x) sum(abs(x)>0)),
                                            nbins = nocc)
                    counts [5,] <- apply(object[,,drop=F],2, function(x) sum(x<0))
                }
                counts [6,] <- apply(object[,,drop=F],2, function(x) sum(abs(x)>0))
                counts [7,] <- apply(object[,,drop=F],2, function(x) length(unique(x[x!=0])))
            }
        }
        if (!is.null(traps)) {
            if (!is.null(usage(traps))) {
                counts[8,] <- apply(usage(traps),2,function(x) sum(x>0))
            }
            else {
                counts[8,] <- rep(ndetector(traps),nocc)
            }
        }

        counts <- as.data.frame(counts)
        dimnames(counts) <- list(c('n','u','f','M(t+1)','losses','detections',
                                   'detectors visited','detectors used'), 1:nocc)
        counts$Total <- apply(counts, 1, sum)
        counts$Total[4] <- counts[4, nocc]
        PSV <-  NULL
        dbar <- NULL

        if (is.null(traps)) {
            trapsum <- NULL
            signalsummary <- NULL
        }
        else {

            if (detector(traps) %in% .localstuff$individualdetectors) {
                if (-diff(counts$Total[1:2]) > 1)
                    PSV <- RPSV(object)
                if (length(dim(object)) == 2)
                    dbar <- dbar(object)
            }
            trapsum <- summary(traps)
            if (detector == 'signal')
                signalsummary <- summary(signal(object))
            if (detector == 'signalnoise')
                signalsummary <- list(signal = summary(signal(object)),
                                      noise = summary(noise(object)),
                                      diffSN = summary(signal(object)-noise(object)))
        }
        zeros <- sum(apply(abs(object),1,sum)==0)
        xyl <- telemetryxy(object)
        if (is.null(xyl))
            telemsummary <- NULL
        else {
            ntelem <- sapply(xyl, nrow)
            nteldet <- if (nrow(object) == 0) 0 else
                sum(apply(abs(object)>0,1,any) [row.names(object) %in% names(xyl)])
            telemsummary <- c(n=length(xyl), ndet=nteldet, min=min(ntelem),
                              max=max(ntelem), mean=mean(ntelem), sd = sd(ntelem))
        }

        ## resightings
        markocc <- markocc(traps(object))
        sighting <- sighting(traps(object))
        if (sighting) {
            sightings <- matrix(0, nrow = 5, ncol = nocc+1, 
                                dimnames = list(c('ID','Not ID','Unmarked',
                                    'Unresolved','Total'), c(1:nocc, 'Total')))
            Tm <- Tm(object)
            Tu <- Tu(object)
            unresolved <- markocc == -1
            sightings[1, c(markocc < 1, FALSE)] <- unlist(counts[6,c(markocc < 1, FALSE)])
            if (!is.null(Tm)) {
                if (length(Tm)==1)
                    sightings[2, 1] <- Tm
                else
                    sightings[2, 1:nocc] <- apply(Tm,2,sum)
            }
            if (!is.null(Tu)) {
                if (is.matrix(Tu))
                    sightings[3, 1:nocc] <- apply(Tu,2,sum)
                else
                    sightings[3, 1] <- Tu
            }                      
            ## 2015-12-15 unresolved = combined 
            if (sum(unresolved)>0){
                sightings[4, unresolved] <- sightings[3,unresolved]    
                sightings[3, unresolved] <- 0    
            }
            sightings[5, 1:nocc] <- apply(sightings[1:4, 1:nocc, drop = FALSE], 2, sum)
                
            sightings[, nocc + 1] <- apply(sightings, 1, sum)
        }
        else sightings <- NULL
        
        temp <- list (
            detector = detector,
            ndetector = nd,
            trapsum = trapsum,
            counts = counts,
            zeros = zeros,
            dbar = dbar,
            RPSV = PSV,
            cutval = cutval,        # signal, signalnoise only
            signalsummary = signalsummary,
            telemsummary = telemsummary,
            sightings = sightings
        )
        class(temp) <- 'summary.capthist'
        temp
    }
}
############################################################################################

counts <- function (CHlist, counts = 'M(t+1)') {
    if (!inherits(CHlist, 'capthist'))
        stop ("require capthist object")
    getc <- function (cnt) {
        getcnt <- function(x, maxocc) {
            temp <- x$counts[cnt,]
            lt <- length(temp)
            matrix(c(temp[-lt], rep(NA, maxocc-lt+1), temp[lt]), nrow = 1)
        }
        if (!is.list(CHlist))
            summary(CHlist)$counts[cnt,]
        else {
            maxocc <- max(sapply(CHlist,ncol))
            abind(lapply(summary(CHlist), getcnt, maxocc), along=1,
                new.names=list(session(CHlist), c(1:maxocc, 'Total')))
        }
    }
    temp <- lapply (counts, getc)
    names(temp) <- counts
    temp
}

print.summary.capthist <- function (x, ...) {
    cat ('Object class     ', 'capthist', '\n')
    nonspatial <- is.null(x$trapsum)
    if (!nonspatial)
        print(x$trapsum, terse=TRUE)
    
    cat ('\nCounts by occasion \n')
    print(x$counts, ...)
    
    if (x$zeros>0)
        cat ('\nEmpty histories : ', x$zeros, '\n')

    if (!nonspatial) {
        if (x$detector %in% c('signal', 'signalnoise')) {
            cat ('Signal threshold ', x$cutval, '\n')
            print (x$signalsummary)
        }
        if (!is.null(x$telemsummary)) {
            cat (x$telemsummary[1], "telemetered animals,", x$telemsummary[2], "detected\n")
            cat (paste(x$telemsummary[3:4], collapse="-"), "locations per animal, mean = ",
                 paste(round(x$telemsummary[5:6],2), collapse=", sd = "), "\n")
        }
    }
    if (!is.null(x$sightings)) {
        cat ('\nSightings by occasion \n')
        print(x$sightings, ...)
        cat ('\n')
    }
    
}
############################################################################################

###############################
## Class : mask
## defines a habitat mask
###############################

subset.mask <- function (x, subset, ...) {

    if (ms(x))
        stop ("subset of multi-session mask not implemented")

    # subset may be numeric index or logical
    temp <- x[subset,,drop=F]
    spacing <- attr(x,'spacing')
    attr(temp,'type')        <- 'subset'
    attr(temp,'meanSD')      <- getMeanSD(temp)
    attr(temp,'area')        <- attr(x, 'area')
    attr(temp,'SLDF')    <- attr(x, 'SLDF')
    attr(temp, 'graph')  <- attr(x, 'graph')
    attr(temp,'spacing')     <- spacing
    if (!is.null(covariates(x))) covariates(temp) <- covariates(x)[subset,,drop=F]
    xl <- range(temp$x) + spacing/2 * c(-1,1)
    yl <- range(temp$y) + spacing/2 * c(-1,1)
    attr(temp,'boundingbox') <- expand.grid(x=xl,y=yl)[c(1,2,4,3),]
    if (!is.null(attr(temp,'OK')))
        attr(temp,'OK') <- attr(temp,'OK')[subset]
    ## 2015-10-18
    ## markingpoints identifies a cutoff: the first mp rows are used (for marking)
    ## this must be reset in a subset operation
    if (!is.null(attr(x,'markingpoints')))
        attr(temp,'markingpoints') <- sum(rep(1,nrow(x))[subset] <= attr(x,'markingpoints'))
    class(temp) <- class(x)
    temp
}
############################################################################################

rbind.mask <- function (...) {
# combine 2 or more mask objects

##    no check for multi-session masks at present
##        stop ('rbind of multi-session mask not implemented')
##
    dropduplicates <- TRUE   ## always
    allargs <- list(...)
    spacing <- attr(allargs[[1]],'spacing')
    area    <- attr(allargs[[1]], 'area')
    check <- function (x) {
        if (!is(x,'mask'))
            stop ("arguments must be mask objects")
        if (attr(x,'spacing') != spacing)
            stop ("arguments must have same 'spacing' attribute")
        if (attr(x,'area') != area)
            stop ("arguments must have same area attribute")
    }
    sapply (allargs, check)
    temp <- rbind.data.frame(...)
    class(temp) <- c('mask', 'data.frame')
    tempcov <- lapply(allargs, covariates)
    covariates(temp) <- do.call(rbind, tempcov)  ## pass list of dataframes

    if (dropduplicates) {
        dupl <- duplicated(temp)
        droppedrows <- sum(dupl)
        if (droppedrows>0) {
            covariates(temp) <- covariates(temp)[!dupl,]
            temp <- temp[!dupl,]
            warning (droppedrows, " duplicate points dropped from mask")
        }
    }

    attr(temp,'type')        <- 'rbind'
    attr(temp,'meanSD')      <- getMeanSD(temp)
    attr(temp,'area')        <- area
    attr(temp,'spacing')     <- spacing
    xl <- range(temp$x) + spacing/2 * c(-1,1)
    yl <- range(temp$y) + spacing/2 * c(-1,1)

    ##  xl <- range(temp[,1])
    ##  yl <- range(temp[,2])

    attr(temp,'boundingbox') <- expand.grid(x=xl,y=yl)[c(1,2,4,3),]
    temp
}
############################################################################################


read.mask <- function (file = NULL, data = NULL, spacing = NULL, columns = NULL, ...)
## 2010-04-11 SS for 'state-space' from SPACECAP
## 2011-05-11 data argument added
## 2011-11-01 re-modelled to add covariates
{
    if (is.null(data) & !is.null(file)) {
        fl <- nchar(file)
        SS <- tolower(substring(file,fl-3,fl)) == '.csv'
        if (SS) {
            data <- read.csv (file)
            if ('HABITAT' %in% names(data))
                data <- data[data$HABITAT == 1,]
        }
        else
            data <- read.table (file, ...)
    }
    else if (is.null(data))
       stop("require one of 'file' or 'data'")
    if (length(dim(data))!=2)
        stop ("require dataframe or matrix for 'data' input to read.mask")

    coln <- colnames(data)
    ixy <- match(c('x','y'), coln)
    if (any(is.na(ixy))) ixy <- 1:2
    mask <- as.data.frame(data[,ixy])
    names(mask) <- c('x', 'y')
    if (any(!apply(mask, 2, is.numeric)))
        stop ("non-numeric x or y coordinates")
    if (any(is.na(mask)))
        stop ("missing value(s) in x or y coordinates")

    class(mask) <- c('mask', 'data.frame')

    ## add covariates
    if (ncol(data) > 2) {
        df <- as.data.frame(data[,-ixy, drop = FALSE])
        if (!is.null(columns)) {
##            if (!all(columns %in% names(mask)))
            ## bug fixed 2014-05-31
            if (!all(columns %in% names(df)))
                stop ("columns missing from input")
            df <- df[,columns, drop=FALSE]
        }
        if (ncol(df)>0)
            covariates(mask) <- df
    }

    if (is.null(spacing))
    {
      sp      <- as.matrix(dist(as.matrix(mask)))
      spacing <- apply(sp,1,function(x) min(x[x>0]))
      spacing <- mean (spacing, na.rm=T)
    }

    area    <- spacing^2 / 10000

    xl <- range(mask$x) + spacing/2 * c(-1,1)
    yl <- range(mask$y) + spacing/2 * c(-1,1)

    attr(mask,'type')    <- 'user'
    attr(mask,'meanSD')  <- getMeanSD(mask)
    attr(mask,'area')    <- area
    attr(mask,'spacing') <- spacing
    attr(mask,'boundingbox') <- expand.grid(x=xl,y=yl)[c(1,2,4,3),]
    attr(mask,'polygon') <- NULL
    mask
}
###############################################################################

summary.mask <- function(object, ...) {

  if (ms(object)) {
      temp <- lapply(object, summary.mask)
      class(temp) <- c('summary.mask', 'list')
      temp
  }
  else {
      if (is.null(object$x) | is.null(object$y))
          stop ("not a valid mask")
      nd <- length(object$x)
      if (length(object$x) != length(object$y))
          stop  ("not a valid mask")

      if (!is.null(covariates(object))) {
          sumcovar <- summary(covariates(object), ...)
      } else sumcovar <- NULL
      if (inherits(object, 'linearmask'))
          maskclass <- 'linearmask'
      else if (inherits(object,'Dsurface'))
          maskclass <- 'Dsurface'
      else if (inherits(object,'Rsurface'))
          maskclass <- 'Rsurface'
      else
          maskclass <- 'mask'

      ## rearranged 2014-09-06
      temp <- list (
        maskclass = maskclass,
        masktype = attr(object, 'type'),
        nmaskpoints = nrow(object),
        xrange = range(object$x),
        yrange = range(object$y),
        meanSD = attr(object, 'meanSD'),
        spacing = attr(object, 'spacing'),
        cellarea = attr(object, 'area'),
        boundingbox = attr(object, 'boundingbox'),
        covar = sumcovar
      )
      class(temp) <- 'summary.mask'
      temp
  }

}
############################################################################################

print.summary.mask <- function (x, ...) {
    if (ms(x)) {
        lapply (x, print.summary.mask)
    }
    else {
      cat ('Object class     ', x$maskclass, '\n')
      cat ('Mask type        ', x$masktype, '\n')
      cat ('Number of points ', x$nmaskpoints, '\n')
      cat ('Spacing m        ', x$spacing, '\n')
      if (is.null(x$cellarea))
          cat ('Total length km  ', x$spacing * x$nmaskpoints / 1000, '\n')
      else {
          cat ('Cell area ha     ', x$cellarea, '\n')
          cat ('Total area ha    ', x$cellarea * x$nmaskpoints, '\n')
      }
      cat ('x-range m        ', x$xrange, '\n')
      cat ('y-range m        ', x$yrange, '\n')
      cat ('Bounding box     ','\n')
      print (x$boundingbox, ...)
      cat ('\n')
      if (!is.null(x$covar)) {
          cat ('Summary of covariates', '\n')
          print(x$covar, ...)
      }
  }
}
############################################################################################

####################################################
## Class : secr
## spatially explicit capture-recapture model fit
####################################################

############################################################################################

trim.secr <- function (object, drop = c('call', 'mask','design','design0'), keep = NULL) {
    trim.default(object, drop = drop, keep = keep)
}
############################################################################################

## 2012-11-14
secrlist <- function(...) {
    dots <- match.call(expand.dots = FALSE)$...
    allargs <- list(...)
    allargs <- lapply(allargs, function(x) if (inherits(x, 'secr')) list(x) else x)
    temp <- do.call(c, allargs)
    ## added 2013-06-06
    if (is.null(names(temp)))
        names(temp) <- paste("secr", 1:length(temp), sep="")
    if (!all(sapply(temp, function(x) inherits(x, 'secr'))))
        stop ("objects must be of class 'secr' or 'secrlist'")
    class(temp) <- 'secrlist'
    temp
}

############################################################################################
## 2014-02-19
## extract method for secrlist objects
## retains class
'[.secrlist' <- function(x, i) {
 y <- NextMethod("[")
 class(y) <- class(x)
 y
}

############################################################################################

coef.secr <- function (object, alpha=0.05, ...) {
    beta   <- object$fit$par
    if (!is.null(object$beta.vcv))
        sebeta <- suppressWarnings(sqrt(diag(object$beta.vcv)))
    else sebeta <- rep(NA, length(beta))
    z <- abs(qnorm(1-alpha/2))
    temp <- data.frame(
        row.names = object$betanames,
        beta    = beta,
        SE.beta = sebeta,
        lcl = beta - z*sebeta,
        ucl = beta + z*sebeta
        )
    attr(temp, 'alpha') <- alpha
    temp
}
############################################################################################

detectpar.default <- function(object, ...) {
    stop ("only for secr models")
}
## byclass option 2013-11-09
detectpar.secr <- function(object, ..., byclass = FALSE) {
    extractpar <- function (temp) {
        if (!is.data.frame(temp))   ## assume list
        {
            if (byclass)
                lapply(temp, extractpar)
            else
                extractpar(temp[[1]])
        }
        else {
            if (!is.data.frame(temp) |
                (nrow(temp) > length(object$link)))
                stop ("unexpected input to detectpar()")

            temp <- temp[, 'estimate', drop = F]
            temp <- split(temp[,1], rownames(temp))
            temp <- c(temp, object$fixed)
            pnames <- parnames(object$detectfn)
            if (object$details$param == 2)
                pnames[1] <- 'esa'
            if (object$details$param %in% c(4,5)) {
                Dindex <- match ('D', names(temp))
                sigmakindex <- match ('sigmak', names(temp))
                cindex <- match ('c', names(temp))
                temp[[sigmakindex]] <- temp[[sigmakindex]] / temp[[Dindex]]^0.5 + temp[[cindex]]
                names(temp)[sigmakindex] <- 'sigma'
            }
            if (object$details$param %in% c(3,5)) {
                a0index <- match ('a0', names(temp))
                sigmaindex <- match ('sigma', names(temp))
                lambda0 <- temp[[a0index]] / 2 / pi / temp[[sigmaindex]]^2 * 10000
                temp[[a0index]] <- if (object$detectfn %in% 14:18) lambda0 else 1-exp(-lambda0)
                names(temp)[a0index] <- if (object$detectfn %in% 0:8) 'g0'
                else if (object$detectfn %in% 14:18) 'lambda0'
                else stop ('invalid combination of param %in% c(3,5) and detectfn')
            }
            temp <- temp[pnames]
            if ((object$detectfn > 9) & (object$detectfn <14))
                temp <- c(temp, list(cutval = object$details$cutval))
            temp
        }
    }
    if (!inherits(object,'secr'))
        stop ("requires 'secr' object")
    temppred <- predict (object, ...)
    if (ms(object)) {
        temp <- lapply(temppred, extractpar)
        names(temp) <- session(object$capthist)
        temp
    }
    else
        extractpar(temppred)
}

############################################################################################

print.secr <- function (x, newdata = NULL, alpha = 0.05, deriv = FALSE, call = TRUE, ...) {

    if (call) {
        cat ('\n')
        print(x$call)
    }

    if (!is.null(x$version)) {
        cat ('secr ', x$version, ', ', x$starttime, '\n', sep='')
    }
    else {   ## for backward compatibility
        cat (x$fitted,'\n')
    }
    cat ('\n')

    print(summary(traps(x$capthist)), terse=TRUE)

    cat ('\n')
    
    ###################
    ## Data description

    if (ms(x$capthist)) {
        print (summary(x$capthist, terse = TRUE))
        Tu <- sapply(x$capthist, function(y) attr(y,'Tu'))
        det <- detector(traps(x$capthist)[[1]])
        xyl <- NULL
    }
    else {
        det <- detector(traps(x$capthist))
        n  <- nrow(x$capthist)     # number caught
        if (length(dim(x$capthist))>2)
            ncapt <- sum(abs(x$capthist))
        else
            ncapt <- sum(abs(x$capthist)>0)

        Tu <- Tu(x$capthist)
        Tm <- Tm(x$capthist)
        unresolvedocc <- markocc(traps(x$capthist)) == -1
        unmarkedocc <- markocc(traps(x$capthist)) == 0
        
        if ('g' %in% x$vars) {
            Groups  <- table(group.factor(x$capthist, x$groups))
            temp <- paste (names(Groups), Groups, collapse=', ', sep='=')
            temp <- paste('(',temp,')', sep='')
        }
        else temp <- ''

        cat ('N animals       : ', n, temp, '\n')
        cat ('N detections    : ', ncapt, '\n')
        if (!is.null(Tu)) {
            cat ('N unmarked sght : ', sum(unlist(Tu)), '\n')
        }    
        if (!is.null(Tm)) 
            cat ('N nonID sghting : ', sum(unlist(Tm)), '\n')
        cat ('N occasions     : ', ncol(x$capthist), '\n')

        xyl <- telemetryxy(x$capthist)
        if (!is.null(xyl)) {
            ## zeros <- sum(apply(abs(object)>0,1,sum)==0)
            ntelem <- sapply(xyl, nrow)
            nteldet <- if (nrow(x$capthist) == 0) 0 else
            sum(apply(abs(x$capthist)>0,1,any) [row.names(x$capthist) %in% names(xyl)])
            ## cat ('Known all-zero  : ', zeros, '\n')
            cat ('Telemetry       : ', length(xyl), 'animals,', nteldet, 'detected\n')
            cat ('Telemetry locns : ', paste(range(ntelem), collapse='-'), 'per animal (mean',
                 round(mean(ntelem),2), 'sd', paste(round(sd(ntelem),2), ')',sep=''), '\n')
        }
        ## cat ('N detectors     : ', ndetector(traps(x$capthist)), '\n')
    }   # end of single-session

    if (det %in% .localstuff$countdetectors) {
        cat ('Count model     :  ')
        if (x$details$binomN == 0) cat ('Poisson \n')
        else if (x$details$binomN == 1) cat ('Binomial, size from usage\n')
        else if (x$details$binomN < 0) cat ('Negative binomial k = ', abs(x$details$binomN), '\n')
        else if (x$details$binomN > 1) cat('Binomial', x$details$binomN, '\n')
    }

    if (!ms(x$capthist)) {
        if (length(maskarea(x$mask))==0)
            cat ('Mask length     : ', masklength(x$mask), 'km \n')
        else
            cat ('Mask area       : ', maskarea(x$mask), 'ha \n')
    }

    ####################
    ## Model description

    ## 2015-03-31
    Npar <- nparameters(x)   ## see utility.R
    AICval <- 2*(x$fit$value + Npar)
    n <- ifelse (ms(x$capthist), sum(sapply(x$capthist, nrow)), nrow(x$capthist))
    AICcval <- ifelse ((n - Npar - 1) > 0,
        2*(x$fit$value + Npar) + 2 * Npar * (Npar+1) / (n - Npar - 1),
        NA)
    cat ('\n')
    cat ('Model           : ', model.string(x$model, x$details$userDfn), '\n')
    if (detector(traps(x$capthist))=='multi') {
        if (!is.null(x$details$param))  ## added 2014-03-13 for compatibility with old secr objects
            if (x$details$param == 1)
                cat ('Gardner, Royle & Wegan parameterisation for multi-catch traps','\n')
    }
    if (!is.null(xyl)) {
        tel12 <- x$details$telemetrytype
        if ((tel12 %in% c('dependent','concurrent')) & x$details$telemetrysigma)
            tel12 <- paste(tel12, 'sigma', sep = ' + ')
        if (tel12 == 'independent')
            tel12 <- paste(tel12, '(sigma only)')
        cat('Telemetry model : ', tel12, '\n')
    }
    ## 2013-06-08
    if (!is.null(x$hcov))
        cat ('Mixture (hcov)  : ', x$hcov, '\n')

    ## 2015-01-16
    if (!is.null(x$details$userdist)) {
        if (is.matrix(x$details$userdist))
            cat ('User distances  :  static (matrix)\n')
        if (is.function(x$details$userdist))
            cat ('User distances  :  dynamic (function)\n')
    }

    cat ('Fixed (real)    : ', fixed.string(x$fixed), '\n')
    cat ('Detection fn    : ', detectionfunctionname(x$detectfn))
    if (!is.null(x$details$normalize))
        if (x$details$normalize)
            cat (', normalized')
    cat ('\n')
    if (!x$CL)
    cat ('Distribution    : ', x$details$distribution, '\n')

    cat ('N parameters    : ', Npar, '\n')
    cat ('Log likelihood  : ', -x$fit$value, '\n')
    cat ('AIC             : ', AICval, '\n')
    cat ('AICc            : ', AICcval, '\n')

    cat ('\n')
    cat ('Beta parameters (coefficients)', '\n')

    print(coef(x), ...)

    if (!is.null(x$fit$hessian)) {
      cat ('\n')
      cat ('Variance-covariance matrix of beta parameters', '\n')
      print (x$beta.vcv, ...)
    }

    # scale newdata covariates... NOT FINISHED 10 05 08
    meanSD <- attr(x$mask,'meanSD')
    if (!is.null(newdata)) {
         for (i in 1:length(newdata)) {
           ind <- match (names(newdata[i]),names(meanSD))
           if (ind>0 & !is.na(meanSD[1,ind]))
             newdata[[i]] <- (newdata[[i]] - meanSD[1,ind]) / meanSD[2,ind]
         }
     }

    cat ('\n')
    cat ('Fitted (real) parameters evaluated at base levels of covariates', '\n')

    if (!is.null(x$realpar))
        print( x$realpar )
    else {
        temp <- predict (x, newdata, type = "response", alpha = alpha)
        nd <- length(temp)
        if (is.data.frame(temp)) print(temp, ...)
        else for (new in 1:nd) {
                cat('\n', names(temp)[new],'\n')
                print(temp[[new]], ...)
            }
    }

    #################################
    # Derived parameters
    #################################
    if (deriv) {

        cat ('\n')
        cat ('Derived parameters', '\n')

        temp <- derived(x, alpha=alpha, se.esa = TRUE)
        nd <- length(temp)
        if (is.data.frame(temp)) print(temp, ...)
        else for (new in 1:nd) {
                cat('\n',names(temp)[new],'\n')
                print(temp[[new]], ...)
            }

    }
    cat ('\n')
}
############################################################################################

oneline.secr <- function (secr) {

    if (ms(secr$capthist)) {
        n <- sum(sapply (secr$capthist, nrow))
        ncapt <- sum(sapply (secr$capthist, function (x) sum(abs(x>0))))
    }
    else  {
        n  <- nrow(secr$capthist)     # number caught
        ncapt <- sum( abs( secr$capthist)>0)
    }

    ## 2015-03-31
    Npar <- nparameters(secr)   ## see utility.R
    AICval <- 2*(secr$fit$value + Npar)
    AICcval <- ifelse ((n - Npar - 1)>0,
        2*(secr$fit$value + Npar) + 2 * Npar * (Npar+1) / (n - Npar - 1),
        NA)

    c (
       model  = model.string(secr$model, secr$details$userDfn),
       detectfn = detectionfunctionname(secr$detectfn),
       npar   = Npar,
       logLik = -secr$fit$value,
       AIC    = round(AICval, 3),
       AICc   = round(AICcval, 3),
       fitted = secr$fitted
    )
}
############################################################################################

logLik.secr <- function(object, ...) {
    npar <- length(object$fit$par)
    structure (-object$fit$value, df = npar, class = 'logLik')
}

############################################################################################

AIC.secr <- function (object, ..., sort = TRUE, k = 2, dmax = 10, criterion = c('AICc','AIC')) {
    allargs <- list(...)
    modelnames <- (c ( as.character(match.call(expand.dots=FALSE)$object),
          as.character(match.call(expand.dots=FALSE)$...) ))
    allargs <- secrlist(object, allargs)
    names(allargs) <- modelnames
    AIC(allargs, sort=sort, k=k, dmax=dmax, criterion=criterion)
}
############################################################################################
############################################################################################

AIC.secrlist <- function (object, ..., sort = TRUE, k = 2, dmax = 10, criterion = c('AICc','AIC')) {

    if (k != 2)
        stop ("AIC.secr defined only for k = 2")

    if (length(list(...)) > 0)
        warning ("... argument ignored in 'AIC.secrlist'")

    if (length(object) > 1) {
        ## check added 2013-10-14
        hcovs <- sapply(object, function(x) if (is.null(x$hcov)) '' else x$hcov)
        if (length(unique(hcovs)) > 1)
            stop ("AIC invalid when models use different hcov")
    }

    criterion <- match.arg(criterion)
    modelnames <- names(object)
    allargs <- object
    if (any(sapply(allargs,class) != 'secr'))
        stop ("components of 'object' must be 'secr' objects")

    output <- data.frame(t(sapply(allargs, oneline.secr)), stringsAsFactors=F)
    for (i in 3:6)
    output[,i] <- as.numeric(output[,i])

    output$delta <- output[,criterion] - min(output[,criterion])
    OK <- abs(output$delta) < abs(dmax)
    sumdelta <- sum(exp(-output$delta[OK]/2))
    output$wt <- ifelse ( OK, round(exp(-output$delta/2) / sumdelta,4), 0)
    row.names(output) <- modelnames
    if (sort) output <- output [order(output[,criterion]),]
    names(output)[7] <- paste('d',criterion,sep='')
    names(output)[8] <- paste(criterion,'wt',sep='')
    if (nrow(output)==1) { output[,8] <- NULL; output[,7] <- NULL}

    output
}
############################################################################################

vcov.secr <- function (object, realnames = NULL, newdata = NULL, byrow = FALSE, ...) {
## return either the beta-parameter variance-covariance matrix
## or vcv each real parameters between points given by newdata (byrow = TRUE)
## or vcv for real parameters at points given by newdata (byrow = TRUE)

    if (is.null(dimnames(object$beta.vcv)))
        dimnames(object$beta.vcv) <- list(object$betanames, object$betanames)

    if (is.null(realnames))
        ## average beta parameters
        return( object$beta.vcv )
    else {
        ## average real parameters
        ## vcv among multiple rows

        if (byrow) {
            ## need delta-method variance of reals given object$beta.vcv & newdata
            if (is.null(newdata))
                newdata <- secr.make.newdata (object)
            nreal <- length(realnames)
            nbeta <- length(object$fit$par)

            rowi <- function (newdatai) {
                reali <- function (beta, rn) {
                    ## real from all beta pars eval at newdata[i,]
                    par.rn <- object$parindx[[rn]]
                    ## 2014-08-19
                    ## mat <- model.matrix(object$model[[rn]], data=newdatai)
                    mat <- general.model.matrix(object$model[[rn]], data=newdatai)
                    lp <- mat %*% matrix(beta[par.rn], ncol = 1)
                    untransform (lp, object$link[[rn]])
                }
                grad <- matrix(nrow = nreal, ncol = nbeta)
                dimnames(grad) <- list(realnames, object$betanames)
                for (rn in realnames)
                    grad[rn,] <- fdHess (pars = object$fit$par, fun = reali, rn = rn)$gradient
                vcv <- grad %*% object$beta.vcv %*% t(grad)
                vcv
            }

            vcvlist <- list(nrow(newdata))
            for (i in 1:nrow(newdata)) vcvlist[[i]] <- rowi(newdata[i,])
            if (length(vcvlist) == 1) vcvlist <- vcvlist[[1]]
            return(vcvlist)
        }
        else {
            newdata <- as.data.frame(newdata)
            rownames <- apply(newdata, 1, function(x) paste(names(newdata), '=', x, sep='',
                                                            collapse=','))
            vcvlist <- list()
            for (rn in realnames) {
                ## temporary fix 2015-09-30
                if (rn == 'pmix')
                    stop("vcov does not work at present when realname == 'pmix'")
                par.rn <- object$parindx[[rn]]
                ## 2014-08-19
                ## mat <- model.matrix(object$model[[rn]], data=newdata)
                mat <- general.model.matrix(object$model[[rn]], data = newdata)
                lp <- mat %*% matrix(object$fit$par[par.rn], ncol = 1)
                real <- untransform (lp, object$link[[rn]])
                real <- as.vector(real)
                ## from Jeff Laake's 'compute.real' in RMark...
                deriv.real <- switch(object$link[[rn]],
                    logit = mat * real * (1-real),
                    log = mat * real,
                    identity = mat,
                    sin = mat * cos(asin(2*real-1))/2)
                vcvlist[[rn]] <- deriv.real %*% object$beta.vcv[par.rn, par.rn] %*% t(deriv.real)
                dimnames(vcvlist[[rn]]) <- list(rownames, rownames)
            }
            names (vcvlist) <- realnames
            return (vcvlist)
        }
        ## DIFFERENT VARIANCE TO secr.lpredictor for sigma because there use se.Xuntransfom
    }
}

############################################################################################

