##############################################################################
## package 'secr'
## fxi.R
## 2014-08-05 : use fxIHP instead of pwuniform
## 2014-08-06 : vector i in fxi.secr, and list output
## 2014-08-27 dist2 optional input to fxIHP set to -1
## 2014-09-03 Ujjwal Kumar's problem: unnecessary minprob in C call
## 2014-09-03 implemented userdist
## 2015-01-31 average all trap locations for default start in fxi.mode
## 2015-02-19 sessnum bug in fxi.contour fixed
## 2015-10-11 mark-resight included
###############################################################################

fxi.secr <- function (object, i = 1, sessnum = 1, X, normal = TRUE, ncores = 1) {

# Return scaled Pr(wi|X).pi(X) for one nominated detection history,
# where X holds coordinates of points
    
    ##--------------------------------------------------------------------
    ## in multi-session case must get session-specific data from lists
    if (ms(object)) {
        session.capthist <- object$capthist[[sessnum]]
        session.traps    <- traps(object$capthist)[[sessnum]]
        session.mask     <- object$mask[[sessnum]]
        session.xy       <- 0
    }
    else {
        session.capthist <- object$capthist
        session.traps    <- traps(object$capthist)
        session.mask     <- object$mask
    }
    ##--------------------------------------------------------------------

    if (missing(X))
        X <- as.matrix(session.mask)
    else
        X <- matrix(unlist(X), ncol = 2)
    beta <- coef(object)$beta
    details <- object$details
    if (is.null(details$param)) details$param <- 0
    Xrealparval  <- makerealparameters (object$design, beta, object$parindx,
                                       object$link, object$fixed)

    ##--------------------------------------------------------------------
    ## validity checks

    if (!all(is.finite(Xrealparval))) {
        cat ('beta vector :', beta, '\n')
        stop ("invalid parameters in 'fxIHP'")
    }

    xylist <- telemetryxy(session.capthist)
    if (!is.null(xylist))
        stop ("fxi.secr is not implemented for telemetry models")

    if (nrow(session.capthist)==0)
        stop ("no data for session ", sessnum)

    if (is.character(i))
        i <- match(i, row.names(session.capthist))  ## convert to numeric index
    if (any(is.na(i) | (i<1) | (i>nrow(session.capthist))))
        stop ("invalid i in fxi.secr")

    if (object$detectfn == 11)
        stop ("fxi.secr does not work with spherical spreading at present")

    if (!is.null(object$groups))
        stop("fxi.secr does not work with groups at present")

    ##--------------------------------------------------------------------
    ## get pdf of location, both pimask and piX
    ## assume uniform for CL = TRUE, or if D homogeneous

    if (is.null(object$model$D))
        D.modelled <- FALSE
    else {
        if (!is.null(object$fixed$D))
            D.modelled <- FALSE
        else
            D.modelled <- (object$model$D != ~1)
    }
    if (D.modelled) {
        predD <- predictDsurface (object)
        if (ms(object))
            predD <- predD[[sessnum]]
        D <- covariates(predD)$D.0  ## does not apply if groups
        pimask <- D / sum(D)   ## vector of probability mass for each mask cell
    }
    else {
        mm <- nrow(session.mask)
        pimask <- rep(1, mm)  ## could be 1/mm, but as we normalise anyway...
    }

    ## fetch predicted density at each new point X
    ## covariates(session.mask) <- data.frame(pi = pimask)
    if (!is.null(covariates(session.mask)))
        covariates(session.mask) <- cbind(data.frame(pi = pimask), covariates(session.mask))
    else
        covariates(session.mask) <- data.frame(pi = pimask)
    ## does this work for linearmask?
    tmpmask <- suppressWarnings(addCovariates(X, session.mask, strict = TRUE))
    piX <- covariates(tmpmask)$pi
    piX[is.na(piX)] <- 0

    #--------------------------------------------------------------------

    nc    <- nrow(session.capthist)
    s     <- ncol(session.capthist)
    m     <- nrow(session.mask)

    dettype <- detectorcode(session.traps)

    if (dettype == 9) {
    # Groups defined by 'animal' covariate
        grp  <- group.factor (session.capthist, 'animal')
        ngrp <- max(1,length(group.levels(session.capthist, 'animal')))
    }
    else {
        grp <- rep(1,nrow(session.capthist))
        ngrp <- 1
    }

    if (dettype %in% c(5,9,12)) {    # signal strength
        session.signal <- signal(session.capthist)
        session.signal <- switch( details$tx,
            log = log(session.signal),
            logit = logit(session.signal),
            identity = session.signal
        )
    }
    else
        session.signal <- 0
    
    
    ##############################################################################
    ## prepare mark-resight data 2015-10-11,15,17
    MRdata <- markresight(session.capthist, session.mask, object$CL, object$fixed,
                          object$details$chat, sessnum) 
    ##############################################################################
    

    ## miscparm is used to package beta parameters that are not modelled
    ## and hence do not have a beta index specified by parindx.
    ## This includes the signal threshold and the mean and sd of noise.

    miscparm <- numeric(4)
    if (object$detectfn %in% c(12,13))       ## experimental signal-noise
        ## fudge: last 2
        miscparm[1:3] <- c(details$cutval,coef(object)[max(unlist(object$parindx))+1:2])
    else if (object$detectfn %in% c(10,11))  ## Dawson & Efford 2009 models
        miscparm[1] <- details$cutval

    if (dettype %in% c(3,6,13)) {            ## polygonX, polygon, telemetry
        k <- table(polyID(session.traps))
        K <- length(k)
        k <- c(k,0)                          ## zero-terminated
        session.xy <- xy(session.capthist)
    }
    else {
        if (dettype %in% c(4,7)) {
            k <- table(transectID(session.traps))
            K <- length(k)
            k <- c(k,0) ## zero terminate
            session.xy <- xy(session.capthist)
        }
        else {
            k <- nrow(session.traps)
            K <- k
            session.xy <- 0
        }
    }
    binomN <- details$binomN
    trps  <- unlist(session.traps, use.names=F)
    usge <- usage(session.traps)
    if (is.null(usge))
        usge <- matrix(1, nrow = K, ncol = s)

    indices <- object$design$PIA[sessnum,1:nc,1:s,1:K,]

    if (is.null(details$userdist)) {
        distmat  <- -1
        distmatX <- -1
    }
    else {
        ## DOES NOT ATTACH COVARIATES FOR noneuc? 2015-02-21
        distmat <- valid.userdist (details$userdist,
                                   detector(session.traps),
                                   xy1 = session.traps,
                                   xy2 = session.mask,
                                   mask = session.mask)
        ## But we also need distances to new points X...
        distmatX <- valid.userdist (details$userdist,
                                   detector(session.traps),
                                   xy1 = session.traps,
                                   xy2 = X,
                                   mask = session.mask)
    }
    fxone <- function (i) {
        temp <- .C('fxIHP', PACKAGE = 'secr',
                   as.integer(i),                 # number of detection history within capthist
                   as.integer(nrow(X)),
                   as.double(X),
                   as.double(piX),
                   as.integer(object$CL),         # 0 = full, 1 = CL
                   as.integer(dettype),           # 0 = multicatch, 1 = proximity, etc
                   as.integer(session.capthist),
                   as.double(unlist(session.xy)),
                   as.double(session.signal),
                   as.integer(grp),
                   as.integer(nc),
                   as.integer(s),
                   as.integer(k),
                   as.integer(m),
                   as.integer(ngrp),
                   as.integer(details$nmix),
                   as.integer(getknownclass(session.capthist, details$nmix, object$hcov)),
                   as.double(trps),
                   as.double(distmat),            ## optional dist2 2014-08-27, 2014-09-03
                   as.double(distmatX),           ## optional dist2 2014-08-27, 2014-09-03
                   as.double(usge),
                   as.integer(MRdata$markocc),
                   as.integer(MRdata$Tu),
                   as.integer(MRdata$Tm),
                   as.double(unlist(session.mask)),
                   as.double(pimask),             ## strictly needed here only if normal>0
                   as.double(Xrealparval),
                   as.integer(nrow(Xrealparval)), # number of rows in lookup table
                   as.integer(indices),           # index of nc,S,K,mix to rows in Xrealparval
                   as.double(miscparm),
                   as.integer(normal),
                   as.integer(object$detectfn),
                   as.integer(binomN),
                   as.double(0),                  ## bugfix 2014-09-03: do not need a floor here
                   value=double(nrow(X)),
                   resultcode=integer(1))
        if (temp$resultcode != 0)
            stop ("error in fxIHP, individual ", i)
        temp$value
    }

    if ((ncores > 1) & (length(i) > 1)) {
        clust <- makeCluster(ncores, methods = FALSE, useXDR = .Platform$endian=='big')
        clusterEvalQ(clust, requireNamespace('secr'))
        output <- parLapply(clust, i, fxone)
        stopCluster(clust)
    }
    else {
        output <- if (length(i) > 1)
            lapply(i, fxone)
        else
            fxone(i)  ## accepts unary operator '-' in fit.mode
    }
    output
}
###############################################################################

fxi2SPDF <- function (x, ID, levels) {
    if (missing(ID))
        ID <- 1:length(x)
    if (missing(levels))
        levels <- names(x[[1]])[names(x[[1]]) != 'mode']
    getxy <- function(one)
        lapply(one[levels], function (xx) Polygon(cbind(xx$x, xx$y)))
    oneanimal <- function (x,id)
        Polygons(x, id)
    xy <- lapply(x[ID], getxy)
    modes <- t(sapply(x[ID], '[[', 'mode') )
    modes <- matrix(unlist(modes), ncol = 2)
    listSrs <- mapply(oneanimal, xy, ID)
    SpP <- SpatialPolygons(listSrs)
    df <- data.frame(modex = modes[,1], modey = modes[,2], row.names = ID)
    SpatialPolygonsDataFrame(SpP, df)
}

##    writeSpatialShape(SPDF, ...)

fxi.contour <- function (object, i = 1, sessnum = 1, border = 100, nx = 64,
    levels = NULL, p = seq(0.1,0.9,0.1), plt = TRUE, add = FALSE, fitmode =
    FALSE, plotmode = FALSE, normal = TRUE, fill = NULL, SPDF = FALSE,  ncores = 1, ...) {
    if (inherits(object$mask, 'linearmask'))
        stop("contouring fxi is not appropriate for linear habitat")
    if (ms(object)) {
        session.traps <- traps(object$capthist)[[sessnum]]
    }
    else {
        session.traps <- traps(object$capthist)
    }
    tempmask <- make.mask (session.traps, border, nx = nx, type = 'traprect')
    xlevels <- unique(tempmask$x)
    ylevels <- unique(tempmask$y)

    fxi <- function (ni) {
        z <- allz[[ni]]
        if (is.null(levels)) {
            temp <- sort(z, decreasing = T)
            cump <- cumsum(temp) / sum(temp)
            levels <- approx (x = cump, y = temp, xout = p)$y
            labels <- p
        }
        else
            labels <- levels

        templines <- contourLines(xlevels, ylevels, matrix(z, nrow = nx), levels = levels)
        ## extra effort to apply correct labels
        getlevels <- function(clines.i) sapply(clines.i, function(q) q$level)
        label.levels <- function (island) {
            which.levels <- match (getlevels(island), levels)
            names(island) <- labels[which.levels]
            island
        }
        templines <- label.levels(templines)

        wh <- which.max(unlist(lapply(templines, function(y) y$level)))
        if (length(templines) > 0) {   ## 2011-04-14
            cc <- templines[[wh]]
            cc <- data.frame(cc[c('x','y')])
            templines$mode <- data.frame(x=mean(cc$x), y=mean(cc$y))

            if (fitmode)
                templines$mode <- fxi.mode(object, sessnum = sessnum,
                                           start = templines$mode, i = ni)
            if (plt) {
                labels <- labels[!is.na(levels)]
                levels <- levels[!is.na(levels)]

                contour (xlevels, ylevels, matrix(z, nrow = nx), add = add,
                         levels = levels, labels = labels, ...)

                ## optional fillin
                if (!is.null(fill)) {
                    z[z < (0.999 * min(levels))] <- NA
                    levels <- rev(c(1,levels,0))
                    .filled.contour(xlevels, ylevels,  matrix(z, nrow = nx), levels= levels,
                                    col = fill)
                }

                if (plotmode) {
                    points(templines$mode, col = 'red', pch = 16)
                }

            }
        }
        templines
    }
    ## sessnum included 2015-02-19 see Ben Augustine email
    allz <- fxi.secr(object, i = i, sessnum = sessnum, X = tempmask, normal = normal,
                     ncores = ncores)
    if (!is.list(allz))
        allz <- list(allz)
    temp <- lapply(1:length(allz), fxi)

    if (SPDF)
        temp <- fxi2SPDF(temp)

    if (plt)
        invisible(temp)
    else
        temp
}

###############################################################################

fxi.mode <- function (object, i = 1, sessnum = 1, start = NULL, ...) {
    if (ms(object))
        session.capthist <- object$capthist[[sessnum]]
    else
        session.capthist <- object$capthist
    start <- unlist(start)
    if (is.null(start)) {
        session.traps <- traps(session.capthist)
        animal <- animalID(session.capthist, names=F) == i

        ## trp <- trap(session.capthist)[animal][1]
        ## start <- unlist(traps(session.capthist)[trp,])
        trp <- trap(session.capthist)[animal]
        start <- sapply(traps(session.capthist)[trp,],mean)
    }
    if (is.character(i))
        i <- match(i, row.names(session.capthist))
    if (is.na(i) | (i<1) | (i>nrow(session.capthist)))
        stop ("invalid i in fxi.secr")

    fn <- function(xy,i) -fxi.secr(object, i = i, sessnum = sessnum, X = xy, normal = FALSE)
    temp <- nlm(f = fn, p = start, i = i, typsize = start, ...)$estimate
    data.frame(x=temp[1], y=temp[2])
}

###############################################################################

## 2014-07-16

## mask if specified should be for a single session
## ... passes newdata df to predict.secr

## upgraded 2015-06-20
## fx.total <- function (object, sessnum = 1, mask = NULL, ncores = 1, ...) {
##     if (ms(object)) stop ("not for ms model")
##     n <- nrow(object$capthist)
##     if (is.null(mask))
##         mask <- if (ms(object)) object$mask[[sessnum]] else object$mask
## 
##     ## sum of individual fxi
##     fxi <- fxi.secr(object, i = 1:n, sessnum = sessnum, X = mask, normal = TRUE, ncores = ncores)
##     fx <- do.call(cbind, fxi)
##     fxt <- apply(fx, 1, sum)
##     fxt <- fxt / getcellsize(mask)   ## length or area 2014-09-10
## 
##     ## predicted uncaught animals
##     D <- predictDsurface(object, mask = mask)
##     if (ms(object)) D <- D[[sessnum]]
##     D <- covariates(D)$D.0
##     CH <- if (ms(object)) object$capthist[[sessnum]] else object$capthist
##     pd <- pdot(X = mask, traps = traps(CH), detectfn = object$detectfn,
##                detectpar = detectpar(object, ...), noccasions = ncol(CH))
##     nct <- D * (1 - pd)
## 
##     covariates(mask) <- data.frame(D.fx = fxt, D.nc = nct, D.sum = fxt + nct)
##     class(mask) <- c('Dsurface', class(mask))
##     mask
## }

fx.total <- function (object, sessnum = 1, mask = NULL, ncores = 1, ...) 
{
    n <- if (ms(object)) nrow(object$capthist[[sessnum]])   
    else nrow(object$capthist)
    if (is.null(mask)) 
        mask <- if (ms(object)) object$mask[[sessnum]]
            else object$mask
    fxi <- fxi.secr(object, i = 1:n, sessnum = sessnum, X = mask, 
        normal = TRUE, ncores = ncores)
    fx <- do.call(cbind, fxi)
    fxt <- apply(fx, 1, sum)
    fxt <- fxt/getcellsize(mask)
    D <- predictDsurface(object, mask = mask)
    D <- covariates(D)$D.0
    CH <- if (ms(object)) 
        object$capthist[[sessnum]]
    else object$capthist
    pd <- pdot(X = mask, traps = traps(CH), detectfn = object$detectfn, 
        detectpar = detectpar(object, ...), noccasions = ncol(CH))
    nct <- D * (1 - pd)
    covariates(mask) <- data.frame(D.fx = fxt, D.nc = nct, D.sum = fxt + 
        nct)
    class(mask) <- c("Dsurface", class(mask))
    mask
}

