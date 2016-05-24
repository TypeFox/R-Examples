#####---------------------------------------------------------------------------
## S3 generic function
getMetric <-
function(x, metric, patID, structure,
     sortBy=c("none", "observed", "patID", "structure", "metric"),
    splitBy=c("none", "patID", "structure", "metric"),
     interp=c("linear", "spline", "ksmooth"), fixed=TRUE, ...) {
    UseMethod("getMetric")
}

## get metrics for one DVH
## patID, structure, sortBy, splitBy are ignored
getMetric.DVHs <-
function(x, metric, patID, structure,
     sortBy=c("none", "observed", "patID", "structure", "metric"),
    splitBy=c("none", "patID", "structure", "metric"),
     interp=c("linear", "spline", "ksmooth"), fixed=TRUE, ...) {

    interp <- match.arg(interp)

    ## D*%  = min radiation of most extremely radiated *% volume
    ## D*cc = min radiation of most extremely radiated *cc volume
    getDose <- function(val, type="relative", unitRef, unitDV) {
        vol <- if(type == "relative") {
            x$dvh[ , "volumeRel"]
        } else if(type == "absolute") {
            x$dvh[ , "volume"]
        }

        ## warn if beyond smallest/largest available volume
        ## like in D0.001cc, D100000cc
        if(isTRUE((val < min(vol)) || (val > max(vol)))) {
            warning(c("Requested ", type, " volume\n", val,
                      " is outside the DVH interval [",
                      formatC(min(vol), format="f", digits=3), ",",
                      formatC(max(vol), format="f", digits=3), "]\n",
                      "NA returned\n"))
            return(NA_real_)
        }

        dose <- if(!is.na(unitDV) && (unitDV == "%")) {   # return relative dose
            x$dvh[ , "doseRel"]
        } else {                         # return absolute dose (default)
            cf <- if(!is.na(unitDV)) {
                getConvFac(paste0(x$doseUnit, "2", unitDV))
            } else {
                1
            }
            cf * x$dvh[ , "dose"]
        }

        if(interp == "linear") {
            ## rule=1 -> no interpolation beyond bounds
            ## ties=max because of DVHs with few unique nodes
            res <- try(approx(vol, dose, val, method="linear", rule=1, ties=max)$y)
            if(!inherits(res, "try-error")) {
                res
            } else {
                NA_real_
            }
        } else if(interp == "ksmooth") {   # kernel smoothing
            bw <- try(KernSmooth::dpill(vol, dose))
            sm <- try(KernSmooth::locpoly(vol, dose, bandwidth=bw,
                                          gridsize=1001L, degree=3))
            if(!inherits(sm, "try-error")) {
                idx <- which.min(abs(sm$x-val))
                sm$y[idx]
            } else {
                NA_real_
            }
        } else if(interp == "spline") {  # does interpolation beyond bounds
            sfun <- try(splinefun(vol, dose, method="monoH.FC", ties=max))
            if(!inherits(sfun, "try-error")) {
                sfun(val)
            } else {
                NA_real_
            }
        }
    }

    ## V*%  = max volume with radiation >= *% of prescribed dose
    ## V*Gy = max volume with radiation >= * Gy / cGy
    getVolume <- function(val, type="relative", unitRef, unitDV) {
        if(type == "relative") {         # given dose is relative
            dose <- x$dvh[ , "doseRel"]

            ## check if max dose is smaller than requested % of prescribed dose
            ## do this here -> while there is approx(..., yright=0)
            ## nothing of that sort exists for spline()
            if(val > 100) {
                warning("max dose is less than requested dose")
                return(0)
            }
        } else if(type == "absolute") {  # given dose is absolute
            cf   <- getConvFac(paste0(x$doseUnit, "2", unitRef))
            dose <- cf * x$dvh[ , "dose"]

            ## check if max dose is smaller than requested % of prescribed dose
            doseMax <- if(!is.na(x$doseMax)) {
                x$doseMax
            } else {
                getDMEAN(x)$doseMax
            }

            if(val > (cf * doseMax)) {
                warning("max dose is less than requested dose")
                return(0)
            }
        }

        vol <- if(!is.na(unitDV) && (unitDV == "CC")) {    # return absolute volume
            x$dvh[ , "volume"]
        } else {                         # return relative volume (default)
            x$dvh[ , "volumeRel"]
        }

        if(interp == "linear") {  # rule=1 -> no interpolation beyond bounds
            res <- try(approx(dose, vol, val, method="linear", rule=1, yright=0)$y)
            if(!inherits(res, "try-error")) {
                res
            } else {
                NA_real_
            }
        } else if(interp == "ksmooth") {   # kernel smoothing
            bw <- try(KernSmooth::dpill(dose, vol))
            sm <- try(KernSmooth::locpoly(dose, vol, bandwidth=bw,
                                          gridsize=1001L, degree=3))
            if(!inherits(sm, "try-error")) {
                idx <- which.min(sm$x-val)
                sm$y[idx]
            } else {
                NA_real_
            }
        } else if(interp == "spline") {  # does interpolation beyond bounds
            sfun <- try(splinefun(dose, vol, val, method="monoH.FC"))
            if(!inherits(sfun, "try-error")) {
                sfun(val)
            } else {
                NA_real_
            }
        }
    }

    ## special metrics that are recognized when prefixed with D, e.g., DMEAN, DEUD
    specMetr <- getSpecialMetrics()

    ## get value for 1 parsed metric
    getVal <- function(pm) {
        if(!pm$valid) {
            return(NA_real_)
        } else if(pm$DV == "D") {              # report a dose
            if(pm$valRef %in% specMetr) {
                cf <- if(!is.na(pm$unitDV)) {
                    getConvFac(paste0(x$doseUnit, "2", pm$unitDV))
                } else {
                    1
                }

                mmmrs <- if(any(c(is.na(c(x$doseMin, x$doseMax, x$doseAvg, x$doseMed, x$doseSD)),
                                  is.null(x$doseMin),
                                  is.null(x$doseMax),
                                  is.null(x$doseAvg),
                                  is.null(x$doseMed),
                                  is.null(x$doseSD))) &&
                            !(pm$valRef %in% c("HI", "EUD", "NTCP", "TCP"))) {
                    getDMEAN(x, interp=interp)
                } else {
                    x
                }

                cf *   if(pm$valRef == "MIN") {
                    if(is.null(x$doseMin) || is.na(x$doseMin)) {
                        mmmrs$doseMin
                    } else { x$doseMin }
                } else if(pm$valRef == "MAX") {
                    if(is.null(x$doseMax) || is.na(x$doseMax)) {
                        mmmrs$doseMax
                    } else { x$doseMax }
                } else if(pm$valRef == "MEAN") {
                    if(is.null(x$doseAvg) || is.na(x$doseAvg)) {
                        mmmrs$doseAvg
                    } else { x$doseAvg }
                } else if(pm$valRef == "MEDIAN") {
                    if(is.null(x$doseMed) || is.na(x$doseMed)) {
                        mmmrs$doseMed
                    } else { x$doseMed }
                } else if(pm$valRef == "SD") {
                    if(is.null(x$doseSD) || is.na(x$doseSD)) {
                        mmmrs$doseSD
                    } else { x$doseSD }
                } else if(pm$valRef == "RX") {
                    if(is.null(x$doseRx) || is.na(x$doseRx)) {
                        NA_real_
                    } else { x$doseRx }
                } else if(pm$valRef == "HI") {
                    getHI(x, ...)$HI
                } else if(pm$valRef == "EUD") {
                    getEUD(x, ...)$EUD
				} else if(pm$valRef == "NTCP") {
					getNTCP(x, ...)$NTCP
				} else if(pm$valRef == "TCP") {
					getTCP(x, ...)$TCP
                } else {
                    warning("Unknown metric reference value")
                    NA_real_
                }
            } else if(pm$unitRef == "%") {
                getDose(pm$valRefNum,   type="relative",
                        unitRef=pm$unitRef, unitDV=pm$unitDV)
            } else if(pm$unitRef == "CC") {
                getDose(pm$valRefNum,   type="absolute",
                        unitRef=pm$unitRef, unitDV=pm$unitDV)
            } else {
                warning("Unknown metric reference value")
                NA_real_
            }
        } else if(pm$DV == "V") {          # report a volume
            if(pm$unitRef == "%") {
                getVolume(pm$valRefNum, type="relative",
                          unitRef=pm$unitRef, unitDV=pm$unitDV)
            } else if(pm$unitRef %in% c("GY", "CGY")) {
                getVolume(pm$valRefNum, type="absolute",
                          unitRef=pm$unitRef, unitDV=pm$unitDV)
            }
        } else {
            warning("Unknown metric")
            NA_real_
        }
    }

    ## parse metric strings into lists of components
    pm  <- parseMetric(metric)
    pmL <- split(pm, f=seq_len(nrow(pm)))
    pmL <- setNames(pmL, vapply(pmL, function(y) y$metric, character(1)))

    ## calculate results for metrics
    Map(getVal, pmL)
}

## getMetric.DVHLst(): return list of separate data frames
## one for each patient, structure, or metric
## x is a list of DVH objects - 1 for each structure or 1 for each patient
getMetric.DVHLst <-
function(x, metric, patID, structure,
     sortBy=c("none", "observed", "patID", "structure", "metric"),
    splitBy=c("none", "patID", "structure", "metric"),
     interp=c("linear", "spline", "ksmooth"), fixed=TRUE, ...) {

    dots <- list(...)

    if(!missing(patID)) {
        ## filter by patID
        p <- trimWS(patID, side="both")
        x <- Filter(function(y) {
            if(fixed) {
                any(y$patID %in% p)
            } else {
                any(grepl(paste(p, collapse="|"), y$patID))
            }
        }, x)

        if(length(x) < 1L) { stop("selected patient not found") }
    }

    if(!missing(structure)) {
        ## filter by structure
        s <- trimWS(structure, side="both")
        x <- Filter(function(y) {
            if(fixed) {
                any(y$structure %in% s)
            } else {
                any(grepl(paste(s, collapse="|"), y$structure))
            }
        }, x)
        if(length(x) < 1L) { stop("selected structure not found") }
    }

    args   <- list(f=getMetric,
                   x=x,
                   metric=list(metric),
                   interp=list(interp))
    res    <- do.call("Map", c(args, dots))
    resDFL <- melt(res, value.name="observed")
    names(resDFL)[names(resDFL) == "L1"] <- "structure"
    names(resDFL)[names(resDFL) == "L2"] <- "metric"

    ## sort if requested
    if(!("none" %in% sortBy)) {
        idx    <- do.call("order", lapply(sortBy, function(y) with(resDFL, get(y))))
        resDFL <- resDFL[idx, ]
    }

    ## split if requested
    if(!("none" %in% splitBy)) {
        splRes <- split(resDFL, lapply(splitBy, function(y) with(resDFL, get(y))))
        ## only keep non-empty components
        resDFL <- Filter(function(y) nrow(y) > 0, splRes)
    }

    ## if result is list of length 1 -> convert to data frame
    if((!is.data.frame(resDFL)) && (length(resDFL) == 1)) {
        resDFL <- resDFL[[1]]
    }

    return(resDFL)
}

## getMetric.DVHLstLst(): return list of separate data frames
## one for each structure, or one for each metric
## x is from many DVH-files -> list of DVH lists
getMetric.DVHLstLst <-
function(x, metric, patID, structure,
     sortBy=c("none", "observed", "patID", "structure", "metric"),
    splitBy=c("none", "patID", "structure", "metric"),
     interp=c("linear", "spline", "ksmooth"), fixed=TRUE, ...) {

    dots <- list(...)

    ## re-organize x into by-patient form if necessary
    isByPat <- attributes(x)$byPat
    if(is.null(isByPat) || !isByPat) {
        x <- reorgByPat(x, byPat=TRUE)
    }

    ## if patIDs are selected, filter those
    if(!missing(patID)) {
        p <- trimWS(patID, side="both")
        x <- Filter(function(y) {
            if(fixed) {
                any(y[[1]]$patID %in% p)
            } else {
                any(grepl(paste(p, collapse="|"), y[[1]]$patID))
            }
        }, x)

        if(length(x) < 1L) { stop("No selected patient found") }
    }

    struct <- if(missing(structure)) {
        NULL
    } else {
        structure
    }

    ## DVH list y is from 1 DVH-file: list of DVHs - many structures for one ID
    collectMetrics <- function(y, metric, structure=NULL, more) {
        ## if structures are selected, filter those
        if(!is.null(structure)) {
            s <- trimWS(structure, side="both")
            y <- Filter(function(z) {
                if(fixed) {
                    any(z$structure %in% s)
                } else {
                    any(grepl(paste(s, collapse="|"), z$structure))
                }
            }, y)

            if(length(y) < 1L) { stop("No selected structure found") }
        }

        args <- list(f=getMetric.DVHs,
                     x=y,
                     metric=list(metric),
                     interp=list(interp))
        do.call("Map", c(args, more))
        #Map(getMetric, y, metric=list(metric), interp=list(interp), dots)
    }

    args <- list(f=collectMetrics,
                 y=x,
                 metric=list(metric),
                 structure=list(struct),
                 more=list(dots))
    res <- do.call("Map", args)
    #res    <- Map(collectMetrics, x, metric=list(metric), structure=list(struct))
    resDFL <- melt(res, value.name="observed")
    names(resDFL)[names(resDFL) == "L1"] <- "patID"
    names(resDFL)[names(resDFL) == "L2"] <- "structure"
    names(resDFL)[names(resDFL) == "L3"] <- "metric"

    ## sort if requested
    if(!("none" %in% sortBy)) {
        idx    <- do.call("order", lapply(sortBy, function(y) with(resDFL, get(y))))
        resDFL <- resDFL[idx, ]
    }

    ## split if requested
    if(!("none" %in% splitBy)) {
        splRes <- split(resDFL, lapply(splitBy, function(y) with(resDFL, get(y))))
        ## only keep non-empty components
        resDFL <- Filter(function(y) nrow(y) > 0, splRes)
    }

    ## if result is list of length 1 -> convert to data frame
    if((!is.data.frame(resDFL)) && (length(resDFL) == 1)) {
        resDFL <- resDFL[[1]]
    }

    return(resDFL)
}
