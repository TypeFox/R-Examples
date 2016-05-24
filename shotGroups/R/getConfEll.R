getConfEll <-
function(xy, level=0.5, dstTarget=100, conversion="m2cm", doRob=TRUE) {
    UseMethod("getConfEll")
}

getConfEll.data.frame <-
function(xy, level=0.5, dstTarget=100, conversion="m2cm", doRob=TRUE) {
    xy <- getXYmat(xy, xyTopLeft=FALSE, relPOA=FALSE)
    NextMethod("getConfEll")
}

getConfEll.default <-
function(xy, level=0.5, dstTarget=100, conversion="m2cm", doRob=TRUE) {
    if(!is.matrix(xy))     { stop("xy must be a matrix") }
    if(!is.numeric(xy))    { stop("xy must be numeric") }
    if(nrow(xy) < 2L)      { stop("we need >= 2 points for confidence ellipse") }
    if(!is.numeric(level)) { stop("level must be numeric") }
    if(level <= 0)         { stop("level must be > 0") }

    ## check if CI level is given in percent
    if(level >= 1) {
        while(level >= 1) { level <- level / 100 }
        warning(c("level must be in (0,1) and was set to ", level))
    }

    ## group center and covariance matrix
    ctr    <- colMeans(xy)               # group center
    covXY  <- cov(xy)                    # covariance matrix (x,y)-coords
    ellRad <- sqrt(eigen(covXY)$values)  # radii error ellipse
    trXY   <- sum(diag(covXY))           # trace of covariance matrix
    detXY  <- det(covXY)                 # determinant

    ## magnification factor to turn error ellipse into confidence ellipse
    N   <- nrow(xy)                      # number of observations
    dfn <- ncol(xy)                      # numerator df
    dfd <- N-1                           # denominator df
    mag <- sqrt(dfn*qf(level, dfn, dfd)) # magnification factor = t-value
    ## radii confidence ellipse
    size <- makeMOA(mag*ellRad, dst=dstTarget, conversion=conversion)
    colnames(size) <- if(dfn == 2L) {
        c("semi-major", "semi-minor")
    } else {
        paste("semi-axis", seq_len(dfn), sep="-")
    }

    ## error ellipse characteristics -> radii = sqrt of eigenvalues
    ## aspect ratio of ellipse = sqrt of kappa condition index
    aspRat <- sqrt(kappa(covXY, exact=TRUE))
    flat   <- 1 - (1/aspRat)             # flattening
    shape  <- c(aspectRatio=aspRat, flattening=flat, trace=trXY, det=detXY)

    haveRob <- if(nrow(xy) < 4L) {
        if(doRob) { warning("We need >= 4 points for robust estimations") }
        FALSE
    } else {
        TRUE
    }                                    # if(nrow(xy) < 4L)

    if(doRob && haveRob) {               # same for robust estimation
        rob       <- robustbase::covMcd(xy)
        ctrRob    <- rob$center
        covXYrob  <- rob$cov
        ellRadRob <- sqrt(eigen(covXYrob)$values)  # radii error ellipse
        trXYrob   <- sum(diag(covXYrob))           # trace of covariance matrix
        detXYrob  <- det(covXYrob)                 # determinant
        ## radii robust confidence ellipse
        sizeRob <- makeMOA(mag*ellRadRob, dst=dstTarget, conversion=conversion)
        colnames(sizeRob) <- colnames(size)

        aspRatRob <- sqrt(kappa(covXYrob, exact=TRUE))   # aspect ratio
        flatRob   <- 1 - (1/aspRatRob)   # flattening
        shapeRob  <- c(aspectRatio=aspRatRob, flattening=flatRob, trace=trXYrob, det=detXYrob)
    } else {
        ## set robust estimates to NULL if not available
        sizeRob  <- NULL
        shapeRob <- NULL
        ctrRob   <- NULL
        covXYrob <- NULL
    }                                    # if(!(doRob & haveRob))

    return(list(ctr=ctr, ctrRob=ctrRob, cov=covXY, covRob=covXYrob,
                size=size, sizeRob=sizeRob,
                shape=shape, shapeRob=shapeRob, magFac=mag))
}
