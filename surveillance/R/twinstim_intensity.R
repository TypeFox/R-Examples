################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Plot the temporal or spatial evolution of the estimated intensity
###
### Copyright (C) 2012-2015 Sebastian Meyer
### $Revision: 1520 $
### $Date: 2015-11-24 15:12:29 +0100 (Die, 24. Nov 2015) $
################################################################################


intensity.twinstim <- function (x, aggregate = c("time", "space"),
    types = 1:nrow(x$qmatrix), tiles, tiles.idcol = NULL)
{
    modelenv <- environment(x)
    
    ## check arguments
    if (is.null(modelenv))
        stop("'x' is missing the model environment\n",
             "  -- re-fit or update() with 'model=TRUE'")
    aggregate <- match.arg(aggregate)
    stopifnot(is.vector(types, mode="numeric"),
              types %in% seq_len(modelenv$nTypes),
              !anyDuplicated(types))

    ## remove (big) x object from current evaluation environment
    qmatrix <- x$qmatrix                # not part of modelenv
    force(types)                        # evaluate types before rm(x)
    rm(x)                               # don't need this anymore

    ##thisenv <- environment()
    ##parent.env(thisenv) <- modelenv     # objects of modelenv become visible
    ## Instead of the above, we do cheap and nasty model unpacking!
    ## safer than the parent.env<- hack (R manual: "extremely dangerous"), and
    ## cleaner than running code inside with(modelenv,...) since assignments
    ## then would take place in modelenv, which would produce garbage
    t0 <- modelenv$t0
    T <- modelenv$T
    histIntervals <- modelenv$histIntervals
    eventTimes <- modelenv$eventTimes
    eventCoords <- modelenv$eventCoords
    eventTypes <- modelenv$eventTypes
    removalTimes <- modelenv$removalTimes
    gridTiles <- modelenv$gridTiles
    gridBlocks <- modelenv$gridBlocks
    ds <- modelenv$ds
    tiaf <- modelenv$tiaf
    tiafpars <- modelenv$tiafpars
    eps.s <- modelenv$eps.s
    siaf <- modelenv$siaf
    siafpars <- modelenv$siafpars
    
    ## endemic component on the spatial or temporal grid
    hInt <- 
        if (modelenv$hash) {
            eta <- drop(modelenv$mmhGrid %*% modelenv$beta)
            if (!is.null(modelenv$offsetGrid)) eta <- modelenv$offsetGrid + eta
            expeta <- exp(unname(eta))
            .beta0 <- rep_len(if (modelenv$nbeta0==0L) 0 else modelenv$beta0,
                              modelenv$nTypes)
            fact <- sum(exp(.beta0[types]))
            if (aggregate == "time") {      # int over W and types by BLOCK
                fact * c(tapply(expeta * modelenv$ds, gridBlocks, sum,
                                simplify = TRUE))
            } else {                        # int over T and types by tile
                fact * c(tapply(expeta * modelenv$dt, gridTiles, sum,
                                simplify = TRUE))
            }
        } else { ## the endemic intensity is 0
            ## but a non-endemic "twinstim" holds no information on 'stgrid':
            ## 'gridBlocks' and 'gridTiles', respectively, are undefined
            NULL
        }

    ## endemic component as a function of time or location
    hIntFUN <- if (modelenv$hash) {
        if (aggregate == "time") {
            function (tp) {
                stopifnot(isScalar(tp))
                if (tp == t0) {
                    hInt[1L]
                } else {
                    starts <- histIntervals$start
                    idx <- match(TRUE, c(starts,T) >= tp) - 1L
                    if (identical(idx, 0L)) { # tp <= t0
                        NA_real_
                    } else { # idx is NA if tp > T
                        block <- histIntervals$BLOCK[idx]
                        hInt[as.character(block)]
                    }
                }
            }
        } else {
            if (!is.null(tiles.idcol)) {
                stopifnot(is(tiles, "SpatialPolygonsDataFrame"))
                row.names(tiles) <- tiles@data[[tiles.idcol]]
            }
            tileLevels <- levels(gridTiles)
            tiles <- check_tiles(tiles, tileLevels,
                                 areas.stgrid = ds[seq_along(tileLevels)],
                                 keep.data = FALSE) # drop data for over-method
            tilesIDs <- row.names(tiles) # = sapply(tiles@polygons, slot, "ID")
            function (xy) {             # works with a whole coordinate matrix
                points <- SpatialPoints(xy, proj4string=tiles@proj4string)
                polygonidxOfPoints <- over(points, tiles)
                tilesOfPoints <- tilesIDs[polygonidxOfPoints]
                hInt[tilesOfPoints]     # index by name
            }
        }
    } else function (...) 0

    ## epidemic component
    eInt <- if (modelenv$hase) {
        qSum_types <- rowSums(qmatrix[,types,drop=FALSE])[eventTypes]
        fact <- qSum_types * modelenv$gammapred
        if (aggregate == "time") {  # as a function of time (int over W & types)
            factS <- fact * modelenv$siafInt
            function (tp) {
                stopifnot(isScalar(tp))
                tdiff <- tp - eventTimes
                infectivity <- qSum_types > 0 & (tdiff > 0) & (removalTimes >= tp)
                if (any(infectivity)) {
                    gsources <- tiaf$g(tdiff[infectivity],
                                       tiafpars,
                                       eventTypes[infectivity])
                    intWj <- factS[infectivity] * gsources
                    sum(intWj)
                } else 0
            }
        } else {  # as a function of location (int over time and types)
            factT <- fact * modelenv$tiafInt
            nEvents <- nrow(eventCoords)
            function (xy) {
                stopifnot(is.vector(xy, mode="numeric"), length(xy) == 2L)
                point <- matrix(xy, nrow=nEvents, ncol=2L, byrow=TRUE)
                sdiff <- point - eventCoords
                proximity <- qSum_types > 0 &
                    .rowSums(sdiff^2, nEvents, 2L) <= eps.s^2
                if (any(proximity)) {
                    fsources <- siaf$f(sdiff[proximity,,drop=FALSE],
                                       siafpars,
                                       eventTypes[proximity])
                    intTj <- factT[proximity] * fsources
                    sum(intTj)
                } else 0
            }
        }
    } else function (...) 0

    ## return component functions
    list(hGrid = hInt, hFUN = hIntFUN, eFUN = eInt,
         aggregate = aggregate, types = types)
}



intensityplot.twinstim <- function (x,
    which = c("epidemic proportion", "endemic proportion", "total intensity"),
    aggregate, types, tiles, tiles.idcol, # arguments of intensity.twinstim;
                                          # defaults are set below
    plot = TRUE, add = FALSE, tgrid = 101, rug.opts = list(),
    sgrid = 128, polygons.args = list(), points.args = list(),
    cex.fun = sqrt, ...)
{
    which <- match.arg(which)
    
    ## set up desired intensities
    cl <- match.call()
    cl <- cl[c(1L, match(names(formals(intensity.twinstim)), names(cl), 0L))]
    cl[[1]] <- as.name("intensity.twinstim")
    components <- eval(cl, envir = parent.frame())
    aggregate <- components$aggregate
    types <- components$types
    
    ## define function to plot
    FUN <- function (tmp) {}
    names(formals(FUN)) <- if (aggregate == "time") "times" else "coords"
    body1 <- if (aggregate == "time") expression(
        hGrid <- sapply(times, components$hFUN, USE.NAMES=FALSE),
        eGrid <- sapply(times, components$eFUN, USE.NAMES=FALSE)
        ) else expression(
            hGrid <- unname(components$hFUN(coords)), # takes whole coord matrix
            eGrid <- apply(coords, 1, components$eFUN)
            )
    body2 <- switch(which,
                    "epidemic proportion" = expression(eGrid / (hGrid + eGrid)),
                    "endemic proportion" = expression(hGrid / (hGrid + eGrid)),
                    "total intensity" = expression(hGrid + eGrid))
    body(FUN) <- as.call(c(as.name("{"), c(body1, body2)))
    
    if (!plot) return(FUN)

    ## plot the FUN
    modelenv <- environment(x)
    dotargs <- list(...)
    nms <- names(dotargs)
    if (aggregate == "time") {
        ## set up grid of x-values (time points where 'which' will be evaluated)
        tgrid <- if (isScalar(tgrid)) {
            seq(modelenv$t0, modelenv$T, length.out=tgrid)
        } else {
            stopifnot(is.vector(tgrid, mode="numeric"))
            sort(tgrid)
        }
        
        ## calculate 'which' on tgrid
        yvals <- FUN(tgrid)
        
        ## plot it
        if(! "xlab" %in% nms) dotargs$xlab <- "time"
        if(! "ylab" %in% nms) dotargs$ylab <- which
        if(! "type" %in% nms) dotargs$type <- "l"
        if(! "ylim" %in% nms) dotargs$ylim <- {
            if (which == "total intensity") c(0,max(yvals)) else c(0,1)
        }
        do.call(if (add) "lines" else "plot", args=c(alist(x=tgrid, y=yvals), dotargs))
        if (is.list(rug.opts)) {
            if (is.null(rug.opts$ticksize)) rug.opts$ticksize <- 0.02
            if (is.null(rug.opts$quiet)) rug.opts$quiet <- TRUE
            eventTimes.types <- modelenv$eventTimes[modelenv$eventTypes %in% types]
            do.call("rug", args = c(alist(x=eventTimes.types), rug.opts))
        }
        invisible(FUN)
    } else {
        tiles <- as(tiles, "SpatialPolygons") # remove potential data for over()
        
        ## set up grid of coordinates where 'which' will be evaluated
        if (isScalar(sgrid)) {
            sgrid <- maptools::Sobj_SpatialGrid(tiles, n = sgrid)$SG
            ## ensure that sgrid has exactly the same proj4string as tiles
            ## since CRS(proj4string(tiles)) might have modified the string
            sgrid@proj4string <- tiles@proj4string
        }
        sgrid <- as(sgrid, "SpatialPixels")
        
        ## only select grid points inside W (tiles)
        sgridTileIdx <- over(sgrid, tiles)
        sgrid <- sgrid[!is.na(sgridTileIdx),]
        
        ## calculate 'which' on sgrid
        yvals <- FUN(coordinates(sgrid))
        sgridy <- SpatialPixelsDataFrame(sgrid, data=data.frame(yvals=yvals),
                                         proj4string=tiles@proj4string)

        ## define sp.layout
        lobjs <- list()
        if (is.list(polygons.args)) {
            nms.polygons <- names(polygons.args)
            if(! "col" %in% nms.polygons) polygons.args$col <- "darkgrey"
            lobjs <- c(lobjs,
                       list(c(list("sp.polygons", tiles, first=FALSE),
                              polygons.args)))
        }
        if (is.list(points.args)) {
            eventCoords.types <- modelenv$eventCoords[modelenv$eventTypes %in% types,,drop=FALSE]
            ## eventCoords as Spatial object with duplicates counted and removed
            eventCoords.types <- SpatialPoints(eventCoords.types,
                                               proj4string = tiles@proj4string,
                                               bbox = tiles@bbox)
            eventCoords.types <- SpatialPointsDataFrame(eventCoords.types,
                data.frame(mult = multiplicity.Spatial(eventCoords.types)))
            eventCoords.types <- eventCoords.types[!duplicated(coordinates(eventCoords.types)),]
            points.args <- modifyList(list(pch=1, cex=0.5), points.args)
            pointcex <- cex.fun(eventCoords.types$mult)
            pointcex <- pointcex * points.args$cex
            points.args$cex <- NULL
            lobjs <- c(lobjs,
                       list(c(list("sp.points", eventCoords.types, first=FALSE,
                                   cex=pointcex), points.args)))
        }
        if ("sp.layout" %in% nms) {
            if (!is.list(dotargs$sp.layout[[1]])) { # let sp.layout be a list of lists
                dotargs$sp.layout <- list(dotargs$sp.layout)
            }
            lobjs <- c(lobjs, dotargs$sp.layout)
            dotargs$sp.layout <- NULL
        }

        ## plotit
        if (add) message("'add'ing is not possible with 'aggregate=\"space\"'")
        if (! "xlim" %in% nms) dotargs$xlim <- bbox(tiles)[1,]
        if (! "ylim" %in% nms) dotargs$ylim <- bbox(tiles)[2,]
        if (! "scales" %in% nms) dotargs$scales <- list(draw = TRUE)
        do.call("spplot", args=c(alist(sgridy, zcol="yvals", sp.layout=lobjs,
                          checkEmptyRC=FALSE), dotargs))
    }
}

## set default arguments for intensityplot.twinstim from intensity.twinstim
formals(intensityplot.twinstim)[names(formals(intensity.twinstim))] <-
    formals(intensity.twinstim)
