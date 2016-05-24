# #############################################################################
# package 'secrlinear'
# methods.R
# 2014-08-28, 2014-09-06, 2014-10-26, 2014-10-30, 2014-11-27
# #############################################################################

plot.linearmask <- function(x, ..., linecol = 'black', label = FALSE,
                            laboffset = c(spacing(x), 0)) {

    defaultarg <- list(col = 'lightgrey', pch = 16, add = FALSE, 
                       lty = 1, lwd = 1, legend = TRUE, pt.cex = 1)
    arg <- list(...)
    arg <- replacedefaults(defaultarg, arg)
    arg$x <- x
     
    linearg <- arg[names(arg) %in% c('lwd','lty')]
    arg <- arg[!names(arg) %in% names(linearg)]
    legend <- arg$legend
    arg$legend <- FALSE   ## suppress plot.mask legend
    
    ## optionally under plot with SLDF
    if (is.null(arg$border)) {
        plot(attr(x, 'SLDF'),  col = linecol, lty = linearg$lty,
             lwd = linearg$lwd, add = arg$add)
        arg$add <- TRUE
    }    
    ## plot.mask
    class(arg$x) <- class(arg$x)[-1]   ## as mask, not linearmask
    legenddata <- do.call(plot, arg)
    
    if (legend & !is.null(arg$covariate)) {
        ncolour <- length(legenddata)
        if (length(arg$col) < ncolour) {
            if (length(arg$col) > 1)
                warning ("too few colours; using terrain.colors(", ncolour, ")")
            arg$col <- terrain.colors(ncolour)   
        }
        legend('right', legend = rev(legenddata), pch = 16, 
               pt.cex = arg$pt.cex, 
               col = rev(arg$col[1:ncolour]),
               title = arg$covariate)
    }    
    
    if (!is.na(linecol))
        plot(attr(x, 'SLDF'),  col = linecol, lty = linearg$lty,
             lwd = linearg$lwd, add = TRUE)
    if (label)
        text (x$x+laboffset[1], x$y+laboffset[2], rownames(x), cex=0.6)
    
    invisible(legenddata)
}
 
# plot.linearmask <- function(x, ..., linecol = 'black', label = FALSE,
#                             laboffset = c(spacing(x), 0), add = FALSE) {
#     arg <- list(...)
#     arg$x <- x
#     if (is.null(arg$col))
#         arg$col <- 'lightgrey'
#     if (is.null(arg$pch))
#         arg$pch <- 16
#     linearg <- arg[names(arg) %in% c('lwd','lty')]
#     if(is.null(linearg$lty))
#         linearg$lty <- 1
#     if(is.null(linearg$lwd))
#         linearg$lwd <- 1
#     arg <- arg[!names(arg) %in% names(linearg)]
#     
#     plot(attr(x, 'SLDF'),  col = linecol, lty = linearg$lty,
#          lwd = linearg$lwd, add = add)
#     do.call(points, arg)
#     plot(attr(x, 'SLDF'),  col = linecol, lty = linearg$lty,
#          lwd = linearg$lwd, add = TRUE)
#     if (label)
#         text (x$x+laboffset[1], x$y+laboffset[2], rownames(x), cex=0.6)
#     
# }

rbind.linearmask <- function (..., cleanskips = TRUE) {

    allargs <- list(...)
    spacing <- attr(allargs[[1]], 'spacing')
    spacingfactor <- attr(allargs[[1]], 'spacingfactor')
    check <- function (x) {
        if (!is(x,'linearmask'))
            stop ("arguments must be linearmask objects")
        if (attr(x,'spacing') != spacing)
            stop ("arguments must have same 'spacing' attribute")
        if (attr(x,'spacingfactor') != spacingfactor)
            stop ("arguments must have same 'spacingfactor' attribute")
    }
    sapply (allargs, check)

    vert <- lapply(allargs, attr, 'SLDF')

    ## must ensure unique ID
    ID <- lapply(vert, function(x) rownames(as(x, 'data.frame')))
    ID2 <- mapply(function(x,i) paste(i,x,sep='.'), ID, 1:length(ID))
    vert <- mapply(spChFIDs, vert, ID2)
    newvert <- do.call(maptools::spRbind, vert)          ## SLDF



    xyl <- lapply(coordinates(newvert), range)
    maskSPDF <- sample.line (newvert, spacing)           ## SPDF
    tmp <- data.frame(maskSPDF)
    mask <- data.frame(coordinates(maskSPDF))            ## dataframe
    names(mask) <- c('x', 'y')

    graph <-  !is.null(attr(allargs[[1]], 'graph'))

    make.linearmask(newvert, spacing, spacingfactor, graph, cleanskips = cleanskips)

}

subset.linearmask <- function (x, subset, LineID, droplinesbeyond = Inf, ...) {

    if (ms(x))
        stop ("subset of multi-session linearmask not implemented")

    # subset may be numeric index or logical
    SLDF <- attr(x, 'SLDF')
    if (!missing(LineID)) {
        if (!missing(subset))
            stop("specify only one of subset and LineID")
        SLDF <- SLDF[LineID, ]
        subset <- covariates(x)$LineID %in% LineID
    }
    temp <- x[subset,,drop=F]
    covariates(temp) <- covariates(x)[subset,,drop=F]

    if (is.finite(droplinesbeyond)) {
        beyond <- function (xy) {
            xy <- do.call(rbind, xy)
            !any(distancetotrap (xy, temp) < droplinesbeyond)
        }
        xyl <- coordinates(SLDF)
        df <- data.frame(SLDF)
        SLDF <- SLDF[!sapply(xyl, beyond), ]
    }
    attr(temp,'type')        <- 'subset'
    attr(temp,'meanSD')      <- getMeanSD(temp)
    attr(temp,'SLDF')        <- SLDF
    attr(temp,'spacing')     <- attr(x,'spacing')
    attr(temp,'spacingfactor')     <- attr(x, 'spacingfactor')

    xyl <- lapply(temp, range)
    names(xyl) <- c('x','y')
    attr(temp,'boundingbox') <- do.call(expand.grid, xyl)[c(1,2,4,3),]

    ## re-form graph if present in first original mask
    if (!is.null(attr(x, 'graph'))) {
        attr(temp, 'graph') <- asgraph(temp)
    }

    if (!is.null(attr(temp,'OK')))
        attr(temp,'OK') <- attr(temp,'OK')[subset]
    class(temp) <- class(x)
    temp
}


plot.linearpopn <- function (x, ..., jitter = 0, plotline = TRUE) {
    mask <- attr(x, 'mask')
    class(x) <- c('popn','data.frame')
    if (!missing(jitter))
        x[,] <- x[,] + jitter * (runif(nrow(x)*2)-0.5) * spacing(mask)
    plot (x, ..., frame = FALSE)
    if (plotline)
        plot (attr(mask,('SLDF')), add = TRUE, col = 'grey')
}
