simRingCount <-
function(xy, target, caliber, unit="cm") {
    UseMethod("simRingCount")
}

simRingCount.data.frame <-
    function(xy, target, caliber, unit="cm") {
        xy <- getXYmat(xy)
    NextMethod("simRingCount")
}

simRingCount.default <-
function(xy, target, caliber, unit="cm") {
    if(!is.matrix(xy))       { stop("xy must be a matrix") }
    if(!is.numeric(xy))      { stop("xy must be numeric") }
    if(ncol(xy) != 2L)       { stop("xy must have two columns") }
    if(!is.numeric(caliber)) { stop("caliber must be numeric") }
    if(caliber <= 0)         { stop("caliber must be > 0") }

    unit <- match.arg(unit, choices=c("cm", "mm", "m", "in", "ft", "yd"))

    ## prepare data
    ## get target definition in requested unit
    trgt <- getTarget(target, unit=unit)

    ## convert caliber to required unit
    convFac <- getConvFac(paste0("mm2", unit))
    calSize <- convFac * caliber

    ## do we need a special simRingCount function?
    if(!is.null(trgt$simRingCount)) {
        fun <- eval(trgt$simRingCount)
        fun(xy, target=trgt, calSize=calSize)
    } else {
        simRingCount_default(xy, target=trgt, calSize=calSize)
    }
}

simRingCount_default <-
function(xy, target, calSize) {
    ## get distance of inner edge of bullet hole to point of aim (0,0)
    ## negative difference -> distance from other side
    dstPOA <- abs(sqrt(rowSums(xy^2)) - calSize/2)

    ## cut with breaks = ring radii
    rings <- if(!is.null(target$countMouche) && target$countMouche) {
        ## with 1st ring (mouche, ring 10 inner sub-division)
        maxAll <- with(target, maxVal+1)
        with(target$inUnit, cut(dstPOA, breaks=c(0, ringR, Inf),
                                labels=c((target$maxVal+1):(target$maxVal-target$nRings+1), 0)),
                                include.lowest=TRUE)
    } else {
        ## except 1st ring (mouche, ring 10 inner sub-division)
        maxAll <- with(target, maxVal)
        with(target$inUnit, cut(dstPOA, breaks=c(0, ringR[-1], Inf),
                                labels=c(target$maxVal:(target$maxVal-target$nRings+1), 0)),
                                include.lowest=TRUE)
    }

    ## convert factor labels to numeric, then sum
    ringCount <- sum(as.numeric(as.character(rings)))  # observed ring count
    ringMax   <- maxAll * nrow(xy)                     # maximum possible

    return(list(count=ringCount, max=ringMax, rings=rings))
}

## TODO
simRingCount_DSUb <-
function(xy, target, calSize) {
    ctrHi <- with(target$inUnit, cbind(0,   ringRV-ringR))
    ctrLo <- with(target$inUnit, cbind(0, -(ringRV-ringR)))

#     ## get convex polygons for all rings
#     getRingPoly <- function(x, ang, i) {
#         pts <- with(x$inUnit, drawDSUOval(c(0, 0),
#             h=ringRV[i]-ringR[i],
#             radius=ringR[i],
#             angle=ang,
#             fg=colsTxt[i],
#             bg=cols[i],
#             plot=FALSE))
#
#         pts[chull(pts), ]
#     }
#
#     ang <- (180/pi)*with(target$inUnit, atan((ringRV-ringR) / ringR))
#     ringPoly <- Map(getRingPoly, list(target), ang, rev(seq_along(target$inUnit$ringR)))

    ## get distance of inner edge of bullet hole to point of aim (0,0)
    ## negative difference -> distance from other side
    dstPOA <- abs(sqrt(rowSums(xy^2)) - calSize/2)

    ## cut with breaks = ring radii
    rings <- if(!is.null(target$countMouche) && target$countMouche) {
        ## with 1st ring (mouche, ring 10 inner sub-division)
        maxAll <- with(target, maxVal+1)
        with(target$inUnit, cut(dstPOA, breaks=c(0, ringR, Inf),
                                labels=c((target$maxVal+1):(target$maxVal-target$nRings+1), 0)),
                                include.lowest=TRUE)
    } else {
        ## except 1st ring (mouche, ring 10 inner sub-division)
        maxAll <- with(target, maxVal)
        with(target$inUnit, cut(dstPOA, breaks=c(0, ringR[-1], Inf),
                                labels=c(target$maxVal:(target$maxVal-target$nRings+1), 0)),
                                include.lowest=TRUE)
    }

    ## convert factor labels to numeric, then sum
    ringCount <- sum(as.numeric(as.character(rings)))  # observed ring count
    ringMax   <- maxAll * nrow(xy)                     # maximum possible

#    return(list(count=ringCount, max=ringMax, rings=rings))
    warning("Simulated ring count not yet available for oval targets")
    return(list(count=NA, max=ringMax, rings=rep(NA_real_, length(rings))))
}

## TODO
simRingCount_DSUa <- simRingCount_DSUb
