###############################################################################
## package 'secr'
## pdot.R
## return net detection probability in 'traps' for home ranges centred at X
## 2010 07 01 alpha detection fn
## 2010 10 09 extended for other detection functions
## 2010 10 10 esa.plot added
## 2010 10 11 esa.plot.secr added
## 2010 10 18 etc. pdot.contour, buffer.contour
## 2010 10 24 tweaking poly
## 2010 11 26 usage
## 2010 11 26 check for ms traps
## 2010 12 19 more careful handling of detectpars
## 2011 01 24 debugged pdotpoly
## 2011 02 06 allow polygonX, transectX
## 2011 06 13 moved spatialscale to utility.R
## 2012 12 24 binomN = 'usage'
## 2014-03-26 pdot.contour and buffer.contour extended to multi-session traps
## 2014-10-17 userdist fixes
## 2014-11-17 more userdist fixes
## 2015-05-15 fill argument for contours
###############################################################################

pdot <- function (X, traps, detectfn = 0, detectpar = list(g0 = 0.2, sigma = 25, z = 1),
                  noccasions = NULL, binomN = NULL, userdist = NULL) {

    ## X should be 2-column dataframe, mask, matrix or similar
    ## with x coord in col 1 and y coord in col 2

    ## added 2010-07-01
    if (is.character(detectfn))
        detectfn <- detectionfunctionnumber(detectfn)
    if ((detectfn > 9) & (detectfn<14) & is.null(detectpar$cutval))
        stop ("requires 'cutval' for detectfn 10:13")
    if (ms(traps))
        stop ("requires single-session traps")

    truncate <- ifelse(is.null(detectpar$truncate), 1e+10, detectpar$truncate)

    detectpars <- unlist(detectpar[parnames(detectfn)])
    if ((detectfn>9) & (detectfn<14))  detectpars <- c(detectpars, detectpar$cutval)

    if (!is.null(usage(traps))) {
        usge <- unlist(usage(traps))
        if (!is.null(noccasions))
            if (noccasions != ncol(usage(traps)))
                warning ("ignoring specified noccasions as it differs from ncol of usage matrix")
        noccasions <- ncol(usage(traps))
    }
    else {
        if (is.null(noccasions))
            stop("must specify noccasions when traps does not have usage attribute")
        usge <- rep(1, ndetector(traps) * noccasions)
    }
    dettype <- detectorcode(traps)
    binomN <- getbinomN (binomN, detector(traps))

    if (!inherits(X, 'mask')) {
        X <- matrix(unlist(X), ncol = 2)
    }
    if (detector(traps) %in% c('polygon','polygonX','transect', 'transectX')) {
        if (!is.null(userdist))
            stop("userdist incompatible with polygon-like detectors")
        k <- table(polyID(traps))   ## also serves transectID
        K <- length(k)              ## number of polygons/transects
        k <-  c(k,0)                ## zero terminate
        temp <- .C('pdotpoly', PACKAGE = 'secr',
            as.double(unlist(X)),
            as.integer(nrow(X)),
            as.double(unlist(traps)),
            as.integer(dettype),
            as.double(usge),
            as.integer(K),
            as.integer(k),
            as.integer(detectfn),   ## hn
            as.double(detectpars),
            as.integer(noccasions),
            as.integer(binomN),
            value = double(nrow(X))
        )
        temp$value
    }
    else {
        #-------------------------------------------------------------
        if (is.null(userdist))
            distmat <- -1
        else {
            distmat <- valid.userdist(userdist,
                                      detector(traps),
                                      xy1 = traps,
                                      xy2 = X,
                                      mask = X)
        }
        #-------------------------------------------------------------
        temp <- .C('pdotpoint', PACKAGE = "secr",
            as.double(unlist(X)),
            as.integer(nrow(X)),
            as.double(unlist(traps)),
            as.double(distmat),                 ## 2014-08-28, 2014-09-01
            as.integer(dettype),
            as.double(usge),
            as.integer(ndetector(traps)),
            as.integer(detectfn),
            as.double(detectpars),
            as.integer(noccasions),
            as.double(truncate^2),
            as.integer(binomN),
            value = double(nrow(X)))
        temp$value
    }
}
############################################################################################

esa.plot <- function (object, max.buffer = NULL, spacing = NULL, max.mask = NULL, detectfn,
                      detectpar, noccasions, binomN = NULL, thin = 0.1, poly = NULL, session = 1,
                      plt = TRUE, as.density = TRUE, n = 1, add = FALSE, overlay = TRUE, ...) {

    if (inherits(object, 'secr')) {
        esa.plot.secr (object, max.buffer, max.mask, thin, poly, session, plt,
                       as.density, add, overlay, ...)
    }
    else {

        if (!inherits(object, 'traps'))
            stop ("requires 'secr' or 'traps' object")
        args <- list(...)
        if(is.null(max.mask)) {
            if (is.null(spacing))
                 spacing <- spacing(object)/3
            max.mask <- make.mask (object, max.buffer, spacing,,, 'trapbuffer', poly)
        }
        detectfn <- valid.detectfn(detectfn)
        binomN <- getbinomN (binomN, detector(object))   ## must now be traps object
        a <- pdot (max.mask, object, detectfn, detectpar, noccasions, binomN)
        d <- distancetotrap(max.mask, object)
        ord <- order(d,a)
        cellsize <-  attr(max.mask, 'spacing')^2/10000
        a <- a[ord]
        output <- data.frame(buffer = d[ord], esa =  cumsum(a) * cellsize,
            density = n /  cumsum(a) / cellsize, pdot = a, pdotmin = cummin(a))
        maxesa <- max(output$esa)
        thinned <- seq(1,  nrow(max.mask), 1/thin)
        output <- output[thinned,]

        if (plt) {
            if (as.density) {
                if (add)
                    lines(output$buffer, n/output$esa, ...)
                else {
                    xlb <- 'Buffer width  m'
                    ylb <- expression(paste('n / esa(buffer)   ',ha^-1))
                    if ('ylim' %in% names(args))
                        plot(output$buffer, n/output$esa, type = 'l',
                            xlab = xlb, ylab = ylb, ...)
                    else  ## clunky!
                        plot(output$buffer, n/output$esa, type = 'l',
                            xlab = xlb, ylab = ylb, ...,
                            ylim= n / maxesa * c(0.9, 1.2))
                }
            }
            else {
                if (add)
                    lines(output$buffer, output$esa, ...)
                else
                    plot(output$buffer, output$esa, type = 'l',
                        xlab = 'Buffer width  m', ylab = 'esa(buffer)  ha', ...)
            }
            invisible(output)
        }
        else output
    }
}

###############################################################################

esa.plot.secr <- function (object, max.buffer = NULL, max.mask = NULL,
    thin = 0.1, poly = NULL, session = 1, plt = TRUE, as.density = TRUE,
    add = FALSE, overlay = TRUE, ...) {

    if (!inherits(object,'secr'))
        stop("require secr object")

    MS <- ms(object)
    if (MS) {
        sessnames <- session(object$capthist)
        ## use alphanumeric session ID
        if (is.numeric(session))
            session <- sessnames[session]
    }

    ## recursive call
    if (MS & (length(session) > 1)) {
        esa.plot.outputs <- vector(mode='list')

        for (i in session) {
            addthisone <- ifelse (add | (overlay & (i != session[1])),
                                  TRUE, FALSE)
            esa.plot.outputs[[i]] <- esa.plot.secr (object, max.buffer,
                max.mask, thin, poly, i, plt, as.density, addthisone,
                overlay, ...)
        }
        if (plt)
            invisible(esa.plot.outputs)
        else
            esa.plot.outputs
    }
    ## not recursive
    else {
        if (MS) {
            ## select one session
            trps <- traps(object$capthist[[session]])
            n <- nrow(object$capthist[[session]])
            nocc <- ncol(object$capthist[[session]])
            spacg <- attr(object$mask[[session]], 'spacing')
            detpar <- detectpar(object)[[session]]
            spscale <- spatialscale(object, object$detectfn, session)
        }
        else {
            trps <- traps(object$capthist)
            n <- nrow(object$capthist)
            nocc <- ncol(object$capthist)
            spacg <- attr(object$mask, 'spacing')
            detpar <- detectpar(object)
            spscale <- spatialscale(object, object$detectfn)
        }
        if (is.null(max.mask)) {
            if (is.null(max.buffer)) {
                if (add)
                    max.buffer <- par()$usr[2]  ## width of existing plot
                else {
                    max.buffer <- 5 * spscale
                }
            }
        }
        binomN <- object$details$binomN
        esa.plot (trps, max.buffer, spacg, max.mask, object$detectfn, detpar,
                  nocc, binomN, thin, poly, session, plt, as.density, n, add, overlay, ...)
    }
}

############################################################################################

pdot.contour <- function (traps, border = NULL, nx = 64, detectfn = 0,
                          detectpar = list(g0 = 0.2, sigma = 25, z = 1),
## noccasions = NULL, binomN = NULL, userdist = NULL,
## no means of passing mask covariates...
                            noccasions = NULL, binomN = NULL,
                          levels = seq(0.1, 0.9, 0.1),
                          poly = NULL, plt = TRUE, add = FALSE, fill = NULL, ...) {
    if (ms(traps)) {
        if (length(noccasions) == 1)
            noccasions <- rep(noccasions,length(traps))
        output <- mapply(pdot.contour, traps, detectpar, noccasions,
                         MoreArgs = list(border = border, nx = nx,
                         detectfn = detectfn, binomN = binomN,
                         levels = levels, poly = poly, plt = plt, add = add, ...))
        if (plt)
            invisible(output)
        else
            output
    }
    else {
        if (is.null(border))
            border <- 5 * spatialscale(detectpar, detectfn)
        tempmask <- make.mask (traps, border, nx = nx, type = 'traprect')
        xlevels <- unique(tempmask$x)
        ylevels <- unique(tempmask$y)
        binomN <- getbinomN (binomN, detector(traps))
        z <- pdot(tempmask, traps, detectfn, detectpar, noccasions, binomN)
        if (!is.null(poly)) {
            OK <- pointsInPolygon(tempmask, poly)
            z[!OK] <- 0
        }
        if (plt) {
            contour (xlevels, ylevels, matrix(z, nrow = nx), add = add, levels = levels, ...)
            
            
            ## optional fillin 2015-05-15
            if (!is.null(fill)) {
                z[z < (0.999 * min(levels))] <- NA
                levels <- c(0,levels,1)
                .filled.contour(xlevels, ylevels,  matrix(z, nrow = nx), levels= levels,
                                col = fill)
            }
            
            
            invisible(contourLines(xlevels, ylevels, matrix(z, nrow = nx), levels = levels))
        }
        else
            contourLines(xlevels, ylevels, matrix(z, nrow = nx), levels = levels)
    }
}
############################################################################################

buffer.contour <- function (traps, buffer, nx = 64, convex = FALSE, ntheta = 100,
                            plt = TRUE, add = FALSE, poly = NULL, fill = NULL, ...) {
    oneconvexbuffer <- function (buffer) {
        temp  <- data.frame(x = apply(expand.grid(traps$x, buffer * cos(theta)),1,sum),
                       y = apply(expand.grid(traps$y, buffer * sin(theta)),1,sum))
        temp <- temp[chull(temp), ]
        temp <- rbind(temp, temp[1,]) ## ensure closed
        rownames(temp) <- NULL
        if (plt)
            lines(temp,...)
        temp
    }
    if (!(inherits(traps, 'traps') | inherits(traps, 'mask')))
        stop ("requires 'traps' or 'mask' object")

    if (ms(traps)) {
        output <- lapply(traps, buffer.contour, buffer = buffer, nx = nx, convex = convex,
               ntheta = ntheta, plt = plt, add = add, poly = poly, ...)
        if (plt)
            invisible(output)
        else
            output
    }
    else {
        if (convex) {
            if (!is.null(poly))
                warning ("'poly' ignored when convex = TRUE")
            ## could use maptools etc. to get intersection?
            theta <- (2*pi) * (1:ntheta) / ntheta
            if (!add & plt)
                plot(traps, border = buffer)
            temp <- lapply(buffer, oneconvexbuffer)
            if (plt)
                invisible (temp)
            else
                temp
        }
        else {
            tempmask <- make.mask (traps, max(buffer)*1.2, nx = nx, type = 'traprect')
            xlevels <- unique(tempmask$x)
            ylevels <- unique(tempmask$y)
            z <- distancetotrap(tempmask, traps)
            if (!is.null(poly)) {
                OK <- pointsInPolygon(tempmask, poly)
                z[!OK] <- 1e20
            }
            if (plt) {
                contour (xlevels, ylevels, matrix(z, nrow = nx), add = add,
                         drawlabels = FALSE, levels = buffer,...)
                if (!is.null(fill)) {
                    z[z < (0.999 * min(levels))] <- NA
                    levels <- c(0,levels,1)
                    .filled.contour(xlevels, ylevels,  matrix(z, nrow = nx), levels= levels,
                                    col = fill)
                }
                invisible(contourLines(xlevels, ylevels, matrix(z, nrow = nx),
                                       levels = buffer))
            }
            else
                contourLines(xlevels, ylevels, matrix(z, nrow = nx),
                             levels = buffer)
        }
    }
}
################################################################################

