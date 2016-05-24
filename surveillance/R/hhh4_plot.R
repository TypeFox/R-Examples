################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Plot-method(s) for fitted hhh4() models
###
### Copyright (C) 2010-2012 Michaela Paul, 2012-2016 Sebastian Meyer
### $Revision: 1687 $
### $Date: 2016-04-01 21:40:25 +0200 (Fre, 01. Apr 2016) $
################################################################################


plot.hhh4 <- function (x,
                       type = c("fitted", "season", "maxEV", "maps", "ri", "neweights"),
                       ...)
{
    stopifnot(x$convergence)
    cl <- sys.call()  # not match.call() because plotHHH4_season() has no 'x'
    ## remove the type argument from the call
    if (is.null(names(cl)) && nargs() > 1L) { # unnamed call plot(x, type)
        cl[[3L]] <- NULL  # remove the second argument
    } else {
        cl$type <- NULL
    }
    cl[[1L]] <- as.name(paste("plotHHH4", match.arg(type), sep="_"))
    eval(cl, envir = parent.frame())
}


###
### Time series of fitted component means and observed counts for selected units
###

plotHHH4_fitted <- function (x, units = 1, names = NULL,
                             col = c("grey85", "blue", "orange"),
                             pch = 19, pt.cex = 0.6, pt.col = 1,
                             par.settings = list(),
                             legend = TRUE, legend.args = list(),
                             legend.observed = FALSE,
                             decompose = NULL, meanHHH = NULL, ...)
{
    if (is.null(units)) units <- seq_len(x$nUnit)
    if (!is.null(names)) stopifnot(length(units) == length(names))
    if (isTRUE(decompose)) decompose <- colnames(x$stsObj)

    ## get decomposed mean
    if (is.null(meanHHH)) {
        meanHHH <- if (is.null(decompose)) {
            meanHHH(x$coefficients, terms.hhh4(x))
        } else {
            decompose.hhh4(x)
        }
    }

    ## check color vector
    col <- if (is.null(decompose) && length(col) == 4) {
        ## compatibility with surveillance < 1.10-0
        pt.col <- col[4L]
        rev(col[-4L])
    } else {
        plotHHH4_fitted_check_col_decompose(col, decompose, dimnames(meanHHH)[[3L]][-1L])
    }

    ## setup graphical parameters
    if (is.list(par.settings)) {
        par.defaults <- list(mfrow = sort(n2mfrow(length(units))),
                             mar = c(4,4,2,0.5)+.1, las = 1)
        par.settings <- modifyList(par.defaults, par.settings)
        opar <- do.call("par", par.settings)
        on.exit(par(opar))
    }

    ## legend options
    if (is.logical(legend)) legend <- which(legend)
    if (!is.list(legend.args)) {
        if (length(legend) > 0)
            warning("ignored 'legend' since 'legend.args' is not a list")
        legend <- integer(0L)
    }
    if (length(legend) > 0) {
        legendidx <- 1L + c(
            if (legend.observed && !is.na(pch)) 0L,
            if (is.null(decompose)) {
                which(c("ne","ar","end") %in% componentsHHH4(x))
            } else seq_along(col))
        default.args <- list(
            x="topright", col=c(pt.col,rev(col))[legendidx], lwd=6,
            lty=c(NA,rep.int(1,length(col)))[legendidx],
            pch=c(pch,rep.int(NA,length(col)))[legendidx],
            pt.cex=pt.cex, pt.lwd=1, bty="n", inset=0.02,
            legend=if (is.null(decompose)) {
                c("observed","spatiotemporal","autoregressive","endemic")[legendidx]
            } else c("observed", rev(decompose), "endemic")[legendidx]
            )
        legend.args <- modifyList(default.args, legend.args)
    }

    ## plot fitted values region by region
    meanHHHunits <- vector(mode="list", length=length(units))
    names(meanHHHunits) <- if (is.character(units)) units else colnames(x$stsObj)[units]
    for(i in seq_along(units)) {
        meanHHHunits[[i]] <- plotHHH4_fitted1(x, units[i], main=names[i],
                                              col=col, pch=pch, pt.cex=pt.cex, pt.col=pt.col,
                                              decompose=decompose, meanHHH=meanHHH, ...)
        if (i %in% legend) do.call("legend", args=legend.args)
    }
    invisible(meanHHHunits)
}

plotHHH4_fitted_check_col_decompose <- function (col, decompose, unitNames)
{
    if (is.null(decompose)) {
        stopifnot(length(col) == 3L)
    } else {
        nUnit <- length(unitNames)
        if (length(col) == nUnit) {
            col <- c("grey85", col)  # first color is for "endemic"
        } else if (length(col) != 1L + nUnit) {
            warning("'col' should be of length ", 1L + nUnit)
            col <- c(col[1L], rep_len(col[-1L], nUnit))
        }
    }
    col
}


### plot estimated component means for a single region

plotHHH4_fitted1 <- function(x, unit=1, main=NULL,
                             col=c("grey85", "blue", "orange"),
                             pch=19, pt.cex=0.6, pt.col=1, border=col,
                             start=x$stsObj@start, end=NULL, xaxis=NULL,
                             xlim=NULL, ylim=NULL, xlab="", ylab="No. infected",
                             hide0s=FALSE, decompose=NULL, meanHHH=NULL)
{
    stsObj <- x$stsObj
    if (is.character(unit) &&
        is.na(unit <- match(.unit <- unit, colnames(stsObj))))
        stop("region '", .unit, "' does not exist")
    if (is.null(main)) main <- colnames(stsObj)[unit]
    if (isTRUE(decompose)) decompose <- colnames(stsObj)

    ## get observed counts
    obs <- observed(stsObj)[,unit]

    ## time range for plotting
    start0 <- yearepoch2point(stsObj@start, stsObj@freq, toleft=TRUE)
    start <- yearepoch2point(start, stsObj@freq)
    tp <- start0 + seq_along(obs)/stsObj@freq # all observation time points
    if (start < start0 || start > tp[length(tp)])
        stop("'start' is not within the time range of 'x$stsObj'")
    end <- if(is.null(end)) tp[length(tp)] else yearepoch2point(end,stsObj@freq)
    stopifnot(start < end)
    tpInRange <- which(tp >= start & tp <= end)            # plot only those
    tpInSubset <- intersect(x$control$subset, tpInRange)   # fitted time points

    ## use time indexes as x-values for use of addFormattedXAxis()
    if (is.list(xaxis)) {
        tp <- seq_along(obs)
        start <- tpInRange[1L]
        end <- tpInRange[length(tpInRange)]
    }

    ## get fitted component means
    meanHHHunit <- if (is.null(decompose)) {
        if (is.null(meanHHH))
            meanHHH <- meanHHH(x$coefficients, terms.hhh4(x))
        sapply(meanHHH, "[", i=TRUE, j=unit)
    } else {
        if (is.null(meanHHH))
            meanHHH <- decompose.hhh4(x)
        if (!setequal(decompose, dimnames(meanHHH)[[3L]][-1L]))
            stop("'decompose' must be (a permutation of) the fitted units")
        meanHHH[,unit,c("endemic",decompose)]
    }
    stopifnot(is.matrix(meanHHHunit), !is.null(colnames(meanHHHunit)),
              nrow(meanHHHunit) == length(x$control$subset))
    meanHHHunit <- meanHHHunit[x$control$subset %in% tpInRange,,drop=FALSE]
    if (any(is.na(meanHHHunit))) { # -> polygon() would be wrong
        ## could be due to wrong x$control$subset wrt the epidemic lags
        ## a workaround is then to set 'start' to a later time point
        stop("predicted mean contains missing values")
    }

    ## check color vector
    col <- if (is.null(decompose) && length(col) == 4L) {
        ## compatibility with surveillance < 1.10-0
        pt.col <- col[4L]
        rev(col[-4L])
    } else {
        plotHHH4_fitted_check_col_decompose(col, decompose, colnames(meanHHHunit)[-1L])
    }
    
    ## establish basic plot window
    if (is.null(ylim)) ylim <- c(0, max(obs[tpInRange],na.rm=TRUE))
    plot(c(start,end), ylim, xlim=xlim, xlab=xlab, ylab=ylab, type="n",
         xaxt = if (is.list(xaxis)) "n" else "s")
    if (is.list(xaxis)) do.call("addFormattedXAxis", c(list(x = stsObj), xaxis))
    title(main=main, line=0.5)

    ## draw polygons
    if (is.null(decompose)) {
        non0 <- which(c("end", "ar", "ne") %in% componentsHHH4(x))
        plotComponentPolygons(
            x = tp[tpInSubset],
            y = meanHHHunit[,c("endemic", "epi.own", "epi.neighbours")[non0],drop=FALSE],
            col = col[non0], border = border[non0], add = TRUE)
    } else {
        non0 <- apply(X = meanHHHunit > 0, MARGIN = 2L, FUN = any)
        plotComponentPolygons(x = tp[tpInSubset], y = meanHHHunit[, non0, drop = FALSE],
                              col = col[non0], border = border[non0], add = TRUE)
    }
    
    ## add observed counts within [start;end]
    ptidx <- if (hide0s) intersect(tpInRange, which(obs > 0)) else tpInRange
    points(tp[ptidx], obs[ptidx], col=pt.col, pch=pch, cex=pt.cex)

    ## invisibly return the fitted component means for the selected region
    invisible(meanHHHunit)
}


### function which does the actual plotting of the polygons

plotComponentPolygons <- function (x, y, col = 1:6, border = col, add = FALSE)
{
    if (!is.vector(x, mode = "numeric") || is.unsorted(x, strictly = TRUE))
        stop("'x' must be a strictly increasing sequence of time points")
    stopifnot(nrow(y <- as.matrix(y)) == (nTime <- length(x)))  # y >= 0
    yc <- if ((nPoly <- ncol(y)) > 1L) {
        apply(X = y, MARGIN = 1L, FUN = cumsum) # nPoly x nTime
    } else t(y)
    
    if (!add) {
        ## establish basic plot window
        plot(range(x), range(yc[nPoly,]), type = "n")
    }
    
    ## recycle graphical parameters
    col <- rep_len(col, nPoly)
    border <- rep_len(border, nPoly)
    
    ## draw polygons
    xpoly <- c(x[1L], x, x[length(x)])
    for (poly in nPoly:1) {
        polygon(x = xpoly, y = c(0, yc[poly, ], 0),
                col = col[poly], border = border[poly])
    }
}


###
### Maps of the fitted mean components averaged over time
###

plotHHH4_maps <- function (x,
    which = c("mean", "endemic", "epi.own", "epi.neighbours"),
    prop = FALSE, main = which, zmax = NULL,
    col.regions = hcl.colors(10),
    labels = FALSE, sp.layout = NULL, ...,
    map = x$stsObj@map, meanHHH = NULL)
{
    which <- match.arg(which, several.ok = TRUE)
    
    ## extract district-specific mean components
    if (is.null(meanHHH)) {
        meanHHH <- meanHHH(x$coefficients, terms.hhh4(x))
    }
    
    ## select relevant components and convert to an array
    meanHHH <- simplify2array(
        meanHHH[c("mean", "endemic", "epi.own", "epi.neighbours")],
        higher = TRUE)
    
    ## convert to proportions
    if (prop) {
        meanHHH[,,-1L] <- meanHHH[,,-1L,drop=FALSE] / c(meanHHH[,,1L])
    }
    
    ## select only 'which' components
    meanHHH <- meanHHH[,,which,drop=FALSE]
    
    ## check map
    map <- as(map, "SpatialPolygonsDataFrame")
    if (!all(dimnames(meanHHH)[[2L]] %in% row.names(map))) {
        stop("'row.names(map)' do not cover all fitted districts")
    }
    
    ## average over time
    comps <- as.data.frame(colMeans(meanHHH, dims = 1))
    
    ## attach to map data
    map@data <- cbind(map@data, comps[row.names(map),,drop=FALSE])

    ## color key range
    if (is.null(zmax)) {
        zmax <- if (prop) {
            ceiling(10*sapply(comps, max))/10
        } else ceiling(sapply(comps, max))
        ## sub-components should have the same color range
        .idxsub <- setdiff(seq_along(zmax), match("mean", names(zmax)))
        zmax[.idxsub] <- suppressWarnings(max(zmax[.idxsub]))
    }

    ## add sp.layout item for district labels
    if (!is.null(layout.labels <- layout.labels(map, labels))) {
        sp.layout <- c(sp.layout, list(layout.labels))
    }
    
    ## produce maps
    grobs <- mapply(
        FUN = function (zcol, main, zmax)
            spplot(map, zcol = zcol, main = main,
                   at = seq(0, zmax, length.out = length(col.regions) + 1L),
                   col.regions = col.regions, sp.layout = sp.layout, ...),
        zcol = names(comps), main = main, zmax = zmax,
        SIMPLIFY = FALSE, USE.NAMES = FALSE)
    if (length(grobs) == 1L) {
        grobs[[1L]]
    } else {
        mfrow <- sort(n2mfrow(length(grobs)))
        gridExtra::grid.arrange(grobs = grobs, nrow = mfrow[1L], ncol = mfrow[2L])
    }
}


###
### Map of estimated random intercepts of a specific component
###

plotHHH4_ri <- function (x, component, labels = FALSE, sp.layout = NULL,
                         gpar.missing = list(col="darkgrey", lty=2, lwd=2),
                         ...)
{
    ranefmatrix <- ranef(x, tomatrix=TRUE)
    if (is.null(ranefmatrix)) stop("model has no random effects")
    stopifnot(length(component) == 1L)
    if (is.na(comp <- pmatch(component, colnames(ranefmatrix))))
        stop("'component' must (partially) match one of ",
             paste(dQuote(colnames(ranefmatrix)), collapse=", "))
    
    map <- as(x$stsObj@map, "SpatialPolygonsDataFrame")
    if (length(map) == 0L) stop("'x$stsObj' has no map")
    map$ranef <- ranefmatrix[,comp][row.names(map)]
    
    if (is.list(gpar.missing) && any(is.na(map$ranef))) {
        sp.layout <- c(sp.layout, 
                       c(list("sp.polygons", map[is.na(map$ranef),]),
                         gpar.missing))
    }
    if (!is.null(layout.labels <- layout.labels(map, labels))) {
        sp.layout <- c(sp.layout, list(layout.labels))
    }
    
    spplot(map[!is.na(map$ranef),], zcol = "ranef",
           sp.layout = sp.layout, ...)
}


###
### Plot the course of the dominant eigenvalue of one or several hhh4-fits
###

plotHHH4_maxEV <- function (...,
                            matplot.args = list(), refline.args = list(),
                            legend.args = list())
{
    objnams <- unlist(lapply(match.call(expand.dots=FALSE)$..., deparse))
    objects <- getHHH4list(..., .names = objnams)

    ## get time points
    epoch <- attr(objects, "epoch")
    start <- attr(objects, "start")
    freq <- attr(objects, "freq")
    start0 <- yearepoch2point(start, freq, toleft=TRUE)
    tp <- start0 + seq_along(epoch) / freq

    ## compute course of dominant eigenvalue for all models
    maxEV <- sapply(objects, getMaxEV, simplify=TRUE, USE.NAMES=TRUE)

    ## line style
    matplot.args <- modifyList(
        list(type="l", col=c(1,2,6,3), lty=c(1,3,2,4), lwd=1.7, cex=1, pch=NULL,
             xlab="", ylab="dominant eigenvalue", ylim=c(0,max(2,maxEV))),
        matplot.args)

    ## main plot
    do.call("matplot", c(list(x=tp, y=maxEV), matplot.args))

    ## add reference line
    if (is.list(refline.args))
        do.call("abline", modifyList(list(h=1, lty=3, col="grey"),
                                     refline.args))

    ## add legend
    if (missing(legend.args) && length(objects) == 1)
        legend.args <- NULL             # omit legend
    if (is.list(legend.args)) {
        legend.args <- modifyList(
            c(list(x="topright", inset=0.02, legend=names(objects), bty="n"),
              matplot.args[c("col", "lwd", "lty", "pch")],
              with(matplot.args, list(pt.cex=cex, text.col=col))),
            legend.args)
        do.call("legend", legend.args)
    }

    ## done
    invisible(maxEV)
}

getMaxEV <- function (x)
{
    Lambda <- createLambda(x)
    if (identical(type <- attr(Lambda, "type"), "zero")) {
        rep.int(0, nrow(x$stsObj))
    } else {
        diagonal <- identical(type, "diagonal")
        vapply(X = seq_len(nrow(x$stsObj)),
               FUN = function (t)
                   maxEV(Lambda(t), symmetric = FALSE, diagonal = diagonal),
               FUN.VALUE = 0, USE.NAMES = FALSE)
    }
}

## generate a function that computes the Lambda_t matrix
createLambda <- function (object)
{
    nTime <- nrow(object$stsObj)
    nUnit <- object$nUnit
    if (identical(componentsHHH4(object), "end")) { # no epidemic components
        zeromat <- matrix(0, nUnit, nUnit)
        Lambda <- function (t) zeromat
        attr(Lambda, "type") <- "zero"
        return(Lambda)
    }
    
    meanHHH <- meanHHH(object$coefficients, terms(object),
                       subset=seq_len(nTime))
    
    W <- getNEweights(object)
    Wt <- if (is.null(W)) {
        NULL
    } else if (is.matrix(W)) {
        function (t) W
    } else {
        function (t) W[,,t]
    }

    type <- NULL
    Lambda <- if (is.null(Wt)) {        # no neighbourhood component
        type <- "diagonal"
        function (t) {
            stopifnot(isScalar(t) && t > 0 && t <= nTime)
            diag(meanHHH$ar.exppred[t,], nUnit, nUnit)
        }
    } else {
        function (t) {
            stopifnot(isScalar(t) && t > 0 && t <= nTime)
            Lambda <- meanHHH$ne.exppred[t,] * Wt(t)
            diag(Lambda) <- diag(Lambda) + meanHHH$ar.exppred[t,]
            Lambda
        }
    }
    attr(Lambda, "type") <- type
    Lambda
}

## determine the dominant eigenvalue of the Lambda matrix
maxEV <- function (Lambda,
                   symmetric = isSymmetric.matrix(Lambda),
                   diagonal = FALSE)
{
    maxEV <- if (diagonal) {
        max(Lambda)  # faster than max(diag(Lambda))
    } else {
        eigen(Lambda, symmetric = symmetric, only.values = TRUE)$values[1L]
    }
    
    ## dominant eigenvalue may be complex
    if (is.complex(maxEV)) {
        if (Im(maxEV) == 0) { # if other eigenvalues are complex
            Re(maxEV)
        } else {
            warning("dominant eigenvalue is complex, using its absolute value")
            abs(maxEV)
        }
    } else {
        maxEV
    }
}



###
### Plot estimated seasonality (sine-cosine terms) of one or several hhh4-fits
### either as multiplicative effect on the 'components' (intercept=FALSE)
### or with intercept=TRUE, which only makes sense if there are no further
### non-centered covariates and offsets.
###

plotHHH4_season <- function (...,
                             components = NULL, intercept = FALSE,
                             xlim = NULL, ylim = NULL,
                             xlab = NULL, ylab = "", main = NULL,
                             par.settings = list(), matplot.args = list(),
                             legend = NULL, legend.args = list(),
                             refline.args = list(), unit = 1)
{
    objnams <- unlist(lapply(match.call(expand.dots=FALSE)$..., deparse))
    objects <- getHHH4list(..., .names = objnams)
    freq <- attr(objects, "freq")
    components <- if (is.null(components)) {
        intersect(c("end", "ar", "ne"), unique(unlist(
            lapply(objects, componentsHHH4), use.names = FALSE)))
    } else {
        match.arg(components, choices = c("ar", "ne", "end", "maxEV"),
                  several.ok = TRUE)
    }

    ## x-axis
    if (is.null(xlim))
        xlim <- c(1,freq)
    if (is.null(xlab))
        xlab <- if (freq==52) "week" else if (freq==12) "month" else "time"

    ## auxiliary function for an argument list "x" with named "defaults" list
    withDefaults <- function(x, defaults)
    {
        if (is.null(x)) defaults else if (is.list(x)) {
            if (is.null(names(x))) {    # x must be complete
                stopifnot(length(x) == length(defaults))
                setNames(x, names(defaults))
            } else modifyList(defaults, x) # x might be a subset of parameters
        } else if (is.atomic(x)) {
            setNames(rep(list(x), length(defaults)), names(defaults))
        } else stop("'", deparse(substitute(x)), "' is not suitably specified")
    }
    
    ## component-specific arguments 
    ylim <- withDefaults(ylim,
                         list(ar=NULL, ne=NULL, end=NULL, maxEV=NULL))
    ylab <- withDefaults(ylab,
                         list(ar=expression(hat(lambda)),
                              ne=expression(hat(phi)),
                              end=expression(hat(nu)),
                              maxEV="dominant eigenvalue"))
    main <- withDefaults(main,
                         list(ar="autoregressive component",
                              ne="spatiotemporal component",    
                              end="endemic component",
                              maxEV="dominant eigenvalue"))
    anyMain <- any(unlist(lapply(main, nchar),
                          recursive=FALSE, use.names=FALSE) > 0)

    ## basic graphical settings
    if (is.list(par.settings)) {
        par.defaults <- list(mfrow=sort(n2mfrow(length(components))),
                             mar=c(4,5,if(anyMain) 2 else 1,1)+.1, las=1)
        par.settings <- modifyList(par.defaults, par.settings)
        opar <- do.call("par", par.settings)
        on.exit(par(opar))
    }

    ## line style
    matplot.args <- modifyList(list(type="l", col=c(1,2,6,3), lty=c(1,3,2,4),
                                    lwd=1.7, cex=1, pch=NULL),
                               matplot.args)

    ## legend options
    if (is.null(legend)) legend <- length(objects) > 1
    if (is.logical(legend)) legend <- which(legend)
    if (!is.list(legend.args)) {
        if (length(legend) > 0)
            warning("ignored 'legend' since 'legend.args' is not a list")
        legend <- integer(0L)
    }
    if (length(legend) > 0) {
        default.args <- c(
            list(x="topright", inset=0.02, legend=names(objects), bty="n"),
            matplot.args[c("col", "lwd", "lty", "pch")],
            with(matplot.args, list(pt.cex=cex, text.col=col))
            )
        legend.args <- modifyList(default.args, legend.args)
    }

    ## plot seasonality in individual model components
    seasons <- list()
    for(comp in setdiff(components, "maxEV")){
        s2 <- lapply(objects, getSeason, component = comp, unit = unit)
        seasons[[comp]] <- exp(vapply(s2, FUN = if (intercept) {
            function (intseas) do.call("+", intseas)
        } else {
            function (intseas) intseas$season  # disregard intercept
        }, FUN.VALUE = numeric(freq), USE.NAMES = TRUE))
        do.call("matplot",              # x defaults to 1:freq
                c(list(seasons[[comp]], xlim=xlim, ylim=ylim[[comp]],
                       xlab=xlab, ylab=ylab[[comp]], main=main[[comp]]),
                  matplot.args))
        if (is.list(refline.args) && !intercept && any(seasons[[comp]] != 1))
            do.call("abline", modifyList(list(h=1, lty=3, col="grey"),
                                         refline.args))
        if (match(comp, components) %in% legend)
            do.call("legend", legend.args)
    }

    ## plot seasonality of dominant eigenvalue
    if ("maxEV" %in% components) {
        seasons[["maxEV"]] <- vapply(objects, FUN = function (obj) {
            getMaxEV_season(obj)$maxEV.season
        }, FUN.VALUE = numeric(freq), USE.NAMES = TRUE)
        do.call("matplot",
                c(list(seasons[["maxEV"]], xlim=xlim,
                       ylim=if (is.null(ylim[["maxEV"]]))
                       c(0,max(2,seasons[["maxEV"]])) else ylim[["maxEV"]],
                       xlab=xlab, ylab=ylab[["maxEV"]],
                       main=main[["maxEV"]]), matplot.args))
        if (is.list(refline.args))
            do.call("abline", modifyList(list(h=1, lty=3, col="grey"),
                                         refline.args))
        if (4 %in% legend) do.call("legend", legend.args)
    }

    ## invisibly return the data that has been plotted
    invisible(seasons)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get estimated intercept and seasonal pattern in the different components
# CAVE: other covariates and offsets are ignored
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getSeason <- function(x, component = c("end", "ar", "ne"), unit = 1)
{
    stopifnot(inherits(x, "hhh4"))
    component <- match.arg(component)
    startseason <- getSeasonStart(x)
    freq <- x$stsObj@freq
    if (is.character(unit)) unit <- match(unit, colnames(x$stsObj))
    
    ## return -Inf is component is not in the model (-> exp(-Inf) = 0)
    if (!component %in% componentsHHH4(x))
        return(list(intercept=-Inf, season=rep.int(-Inf, freq)))

    ## get the intercept
    est <- fixef(x, reparamPsi=FALSE)
    intercept <- unname(est[grep(paste0("^", component, "\\.(1|ri)"), names(est))])
    if (length(intercept) == 0) {
        intercept <- 0 # no intercept (not standard)
    } else if (length(intercept) > 1) { # unit-specific intercepts
        if (length(intercept) != ncol(x$stsObj))
            stop(component,"-component has incomplete unit-specific intercepts")
        intercept <- intercept[unit]
        if (is.na(intercept)) stop("the specified 'unit' does not exist")
    }

    ## get seasonality terms (relying on sin(2*pi*t/52)-kind coefficient names)
    coefSinCos <- est[grep(paste0("^",component, "\\.(sin|cos)\\("), names(est))]
    if (unitspecific <- length(grep(").", names(coefSinCos), fixed=TRUE))) {
        if (unitspecific < length(coefSinCos))
            stop("cannot handle partially unit-specific seasonality")
        coefSinCos <- coefSinCos[grep(paste0(").",colnames(x$stsObj)[unit]),
                                      names(coefSinCos), fixed=TRUE)]
        ## drop .unitname-suffix since non-syntactic (cannot reformulate())
        names(coefSinCos) <- sub("\\)\\..+$", ")", names(coefSinCos))
    }
    if (length(coefSinCos)==0)
        return(list(intercept=intercept, season=rep.int(0,freq)))
    fSinCos <- reformulate(
        sub(paste0("^",component,"\\."), "", names(coefSinCos)),
        intercept=FALSE)
    mmSinCos <- model.matrix(fSinCos,
                             data=data.frame(t=startseason-1 + seq_len(freq)))

    ## Done
    list(intercept=intercept, season=as.vector(mmSinCos %*% coefSinCos))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# compute dominant eigenvalue of Lambda_t
# CAVE: no support for Lambda_it
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getMaxEV_season <- function (x)
{
    stopifnot(inherits(x, "hhh4"))
    nUnits <- x$nUnit
    freq <- x$stsObj@freq
    components <- componentsHHH4(x)

    ## CAVE: this function ignores epidemic covariates/offsets
    ##       and unit-specific seasonality
    if (nUnits > 1L && any(c("ar", "ne") %in% components)) {
        compOK <- vapply(x$control[c("ar","ne")], FUN = function (comp) {
            terms <- terms(x)$terms
            epiterms <- terms[,terms["offsetComp",] %in% seq_len(2L),drop=FALSE]
            identical(as.numeric(comp$offset), 1) &&
                length(all.vars(removeTimeFromFormula(comp$f))) == 0L &&
                    all(!unlist(epiterms["unitSpecific",]))
        }, FUN.VALUE = TRUE, USE.NAMES = FALSE)
        if (any(!compOK))
            warning("epidemic components have (unit-specific) ",
                    "covariates/offsets not accounted for;\n",
                    "  use getMaxEV() or plotHHH4_maxEV()")
    }
    
    ## global intercepts and seasonality
    s2.lambda <- getSeason(x, "ar")
    s2.phi <- getSeason(x, "ne")
    
    ## unit-specific intercepts
    ris <- ranef(x, tomatrix=TRUE)
    ri.lambda <- ris[,pmatch("ar.ri", colnames(ris), nomatch=0L),drop=TRUE]
    if (length(ri.lambda) == 0L) ri.lambda <- rep.int(0, nUnits)
    ri.phi <- ris[,pmatch("ne.ri", colnames(ris), nomatch=0L),drop=TRUE]
    if (length(ri.phi) == 0L) ri.phi <- rep.int(0, nUnits)

    ## get neighbourhood weights as a function of time
    W <- getNEweights(x)  # NULL, matrix or 3-dim array
    if (!is.null(W) && !is.matrix(W))
        stop("neighbourhood weights are time-varying; ",
            # and thus probably changing within or across seasons
             "use getMaxEV() or plotHHH4_maxEV()")

    ## create the Lambda_t matrix
    createLambda <- function (t)
    {
        Lambda <- if ("ne" %in% components) {
            exp(s2.phi$intercept + ri.phi + if(t==0) 0 else s2.phi$season[t]) * W
        } else matrix(0, nUnits, nUnits)
        if ("ar" %in% components) {
            diag(Lambda) <- diag(Lambda) + exp(s2.lambda$intercept + ri.lambda +
                if(t==0) 0 else s2.lambda$season[t])
        }
        Lambda
    }

    ## do this for t in 0:freq
    diagonal <- !("ne" %in% components)
    .maxEV <- function (t) {
        maxEV(createLambda(t), symmetric = FALSE, diagonal = diagonal)
    }
    maxEV.const <- .maxEV(0)
    maxEV.season <- if (all(c(s2.phi$season, s2.lambda$season) %in% c(-Inf, 0))) {
        rep.int(maxEV.const, freq)
    } else {
        vapply(X = seq_len(freq), FUN = .maxEV, FUN.VALUE = 0, USE.NAMES = FALSE)
    }

    ## Done
    list(maxEV.season = maxEV.season,
         maxEV.const  = maxEV.const,
         Lambda.const = createLambda(0))
}


## Determine the time point t of the start of a season in a hhh4() fit.
## If \code{object$stsObj@start[2] == 1}, it simply equals
## \code{object$control$data$t[1]}. Otherwise, the \code{stsObj} time series
## starts within a year (at sample \code{s}, say) and the beginning of
## the next season is
## \code{object$control$data$t[1] + object$stsObj@freq - s + 1}.
getSeasonStart <- function (object)
{
    if ((startsample <- object$stsObj@start[2]) == 1) {
        object$control$data$t[1L]
    } else {
        object$control$data$t[1L] + object$stsObj@freq-startsample + 1
    }
}



###
### plot neighbourhood weight as a function of distance (neighbourhood order)
###

plotHHH4_neweights <- function (x, plotter = boxplot, ...,
                                exclude = 0, maxlag = Inf)
{
    plotter <- match.fun(plotter)
    
    ## orders of neighbourhood (o_ji)
    nbmat <- neighbourhood(x$stsObj)
    if (all(nbmat %in% 0:1)) {
        message("'neighbourhood(x$stsObj)' is binary; ",
                "computing neighbourhood orders ...")
        nbmat <- nbOrder(nbmat, maxlag=maxlag)
    }

    ## extract (estimated) weight matrix (w_ji)
    W <- getNEweights(x)
    if (is.null(W)) {  # if no spatio-temporal component in the model
        W <- nbmat
        W[] <- 0
    }

    ## draw the boxplot
    Distance <- factor(nbmat, exclude = exclude)
    notexcluded <- which(!is.na(Distance))
    Distance <- Distance[notexcluded]
    Weight <- W[notexcluded]
    plotter(Weight ~ Distance, ...)
}



###
### auxiliary functions
###

yearepoch2point <- function (yearepoch, frequency, toleft=FALSE)
    yearepoch[1L] + (yearepoch[2L] - toleft) / frequency


getHHH4list <- function (..., .names = NA_character_)
{
    objects <- list(...)
    if (length(objects) == 1L && is.list(objects[[1L]]) &&
        inherits(objects[[1L]][[1L]], "hhh4")) {
        ## ... is a single list of fits
        objects <- objects[[1L]]
        if (is.null(names(objects))) names(objects) <- seq_along(objects)
    } else {
        names(objects) <- if (is.null(names(objects))) .names else {
            ifelse(nzchar(names(objects)), names(objects), .names)
        }
    }
    if (!all(sapply(objects, inherits, what="hhh4")))
        stop("'...' must consist of hhh4()-fits only")

    ## check common epoch, start and frequency and append them as attributes
    epoch <- unique(t(sapply(objects, function(x) x$stsObj@epoch)))
    if (nrow(epoch) > 1)
        stop("supplied hhh4-models obey different 'epoch's")
    attr(objects, "epoch") <- drop(epoch)
    start <- unique(t(sapply(objects, function(x) x$stsObj@start)))
    if (nrow(start) > 1)
        stop("supplied hhh4-models obey different start times")
    attr(objects, "start") <- drop(start)
    freq <- unique(sapply(objects, function(x) x$stsObj@freq))
    if (length(freq)>1)
        stop("supplied hhh4-models obey different frequencies")
    attr(objects, "freq") <- freq
    
    ## done
    return(objects)
}
