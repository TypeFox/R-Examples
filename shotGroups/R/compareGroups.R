compareGroups <-
function(DF, plots=TRUE, xyTopLeft=TRUE, ABalt=c("two.sided", "less", "greater"),
         Walt=c("two.sided", "less", "greater"), CEPtype="CorrNormal",
         CEPlevel=0.5, CIlevel=0.95, conversion="m2cm") {
    if(!is.data.frame(DF))    { stop("DF must be a data frame") }
    if(!is.numeric(CEPlevel)) { stop("CEPlevel must be numeric") }
    if(CEPlevel <= 0)         { stop("CEPlevel must be > 0") }
    if(!is.numeric(CIlevel))  { stop("CIlevel must be numeric") }
    if(CIlevel <= 0)          { stop("CIlevel must be > 0") }

    CEPtype <- match.arg(CEPtype,
                         choices=c("CorrNormal", "GrubbsPearson", "GrubbsLiu",
                                   "GrubbsPatnaik", "Rayleigh", "Krempasky",
                                   "Ignani", "RMSE", "Ethridge", "RAND", "Valstar"), several.ok=FALSE)

    ## check if CEP / CI level is given in percent
    if(CEPlevel >= 1) {
        while(CEPlevel >= 1) { CEPlevel <- CEPlevel / 100 }
        warning(c("CEPlevel must be in (0,1) and was set to ", CEPlevel))
    }

    if(CIlevel >= 1) {
        while(CIlevel >= 1) { CIlevel <- CIlevel / 100 }
        warning(c("CIlevel must be in (0,1) and was set to ", CIlevel))
    }

    ## convert DF names to lower case
    DF <- setNames(DF, tolower(names(DF)))

    #####-----------------------------------------------------------------------
    ## make sure DF has the required variable names and at least two groups
    varNames <- names(DF)                # what variables are present
    needsSer <- "series"                 # required
    needsXY1 <- c("point.x", "point.y")  # coordinates must have this name
    needsXY2 <- c("x", "y")              # or this
    wantsDst <- "distance"               # useful
    wantsAIM <- c("aim.x", "aim.y")      # useful
    hasSer   <- needsSer %in% varNames   # required ones we have
    hasXY1   <- needsXY1 %in% varNames   # coordinates we have
    hasXY2   <- needsXY2 %in% varNames
    hasDst   <- wantsDst %in% varNames   # useful ones we have
    hasAIM   <- wantsAIM %in% varNames   # useful ones we have

    if(!all(hasSer)) {
        stop(c("The data frame is missing variable\n",
               paste(needsSer[!hasSer], collapse=" ")))
    }

    if(!xor(all(hasXY1), all(hasXY2))) {
        stop("Coordinates must be named either X, Y or Point.X, Point.Y")
    }

    if(!all(hasDst)) {
        warning(c("The data frame is missing variable\n",
                  paste(wantsDst[!hasDst], collapse=" "), "\n",
                  "Distance is assumed to be 100"))
        DF$distance <- 100
    }

    if(!all(hasAIM)) {
        warning(c("The data frame is missing variable(s)\n",
                  paste(wantsAIM[!hasAIM], collapse=" "), "\n",
                  "Point of Aim is assumed to be in (0,0)"))
        DF$aim.x <- 0
        DF$aim.y <- 0
    }

    #####-----------------------------------------------------------------------
    ## prepare data
    res <- vector("list", 0)               # empty list to later collect the results
    if(!is.factor(DF$series)) {            # make sure series is a factor
        DF$series <- as.factor(DF$series)
    } else {
        DF$series <- droplevels(DF$series) # remove all non-used factor levels
    }

    ## check if we have enough groups and points per group
    if(nlevels(DF$series) < 2L)            { stop("We need >= 2 groups for a comparison") }
    if(any(xtabs(~ series, data=DF) < 2L)) { stop("We need >= 2 points in each group") }

    ## prepare data: get (x,y)-coords relative to point of aim as matrix
    xy <- getXYmat(DF, xyTopLeft=xyTopLeft)
    DF <- cbind(DF, xy)
    dstTarget <- tapply(DF$distance, DF$series, mean)  # distances to target

    ## for each group extract the old (x,y)-coords as a matrix
    extractXY <- function(x) {
        DFsub <- x[ , which(!(names(x) %in% c("x", "y")))]
        getXYmat(DFsub, xyTopLeft=xyTopLeft)
    }

    xyL  <- lapply(split(DF, DF$series), extractXY)
    nS   <- length(xyL)                   # total number of series
    nObs <- vapply(xyL, nrow, integer(1)) # number of obs per series
    names(xyL) <- levels(DF$series)

    ## to determine axis limits later, collect all results in a vector
    axisLimsX <- numeric(0)
    axisLimsY <- numeric(0)

    #####-----------------------------------------------------------------------
    ## location measures
    res$ctr <- sapply(xyL, colMeans)     # group centers
    distPOA <- sqrt(colSums(res$ctr^2))  # distances to point of aim
    distPOAmoa  <- Map(makeMOA, distPOA, dst=dstTarget, conversion=conversion)
    res$distPOA <- do.call("cbind", distPOAmoa)

    ## multivariate location test for equal group centers (relative to POA)
    res$MANOVA <- anova(lm(cbind(x, y) ~ series, data=DF), test="Wilks")

    #####-----------------------------------------------------------------------
    ## shape measures
    ## correlation matrices for x- and y-coords
    res$corXY <- lapply(xyL, cor)

    #####-----------------------------------------------------------------------
    ## spread measures
    ## standard deviations for x- and y-coords with parametric CIs
    alpha <- 1-CIlevel
    getSDxyCI <- function(x) {
        N    <- nrow(x)
        sdXY <- sqrt(diag(cov(x)))       # standard deviations

        ## and their parametric confidence intervals
        sdXci  <- sqrt((N-1)*sdXY["x"]^2 / qchisq(c(1-(alpha/2), alpha/2), N-1))
        sdYci  <- sqrt((N-1)*sdXY["y"]^2 / qchisq(c(1-(alpha/2), alpha/2), N-1))
        sdXYci <- c(sdXci, sdYci)

        setNames(sdXYci, c("sdX (", "sdX )", "sdY (", "sdY )"))
    }

    ## sd and sd CIs as separate lists for unit of measurement MOA, SMOA, mrad
    sdXY       <- lapply(xyL, function(x) sqrt(diag(cov(x)))) # standard deviations
    sdXYci     <- lapply(xyL, getSDxyCI) # confidence intervals
    res$sdXY   <- Map(makeMOA, sdXY,   dst=dstTarget, conversion=conversion)
    res$sdXYci <- Map(makeMOA, sdXYci, dst=dstTarget, conversion=conversion)

    ## (mean) distances to group center
    dstCtrL    <- lapply(xyL, getDistToCtr) # distances to group center
    meanDstCtr <- lapply(dstCtrL, mean)
    meanDstCtrMOA <- Map(makeMOA, meanDstCtr, dst=dstTarget, conversion=conversion)
    res$meanDistToCtr <- do.call("cbind", meanDstCtrMOA)

    ## maximum pairwise distance = maximum group spread
    maxPD      <- lapply(xyL, getMaxPairDist)   # max pairwise distance
    maxSpread  <- lapply(maxPD, function(x) { x$d } )
    maxPDidx   <- sapply(maxPD, function(x) { x$idx } )
    maxSpreadL <- Map(makeMOA, maxSpread, dst=dstTarget, conversion=conversion)
    res$maxPairDist <- do.call("cbind", maxSpreadL)

    ## bounding box figure of merit and diagonal
    ## bbs     <- lapply(xyL, getBoundingBox)   # bounding boxes
    bbs     <- lapply(xyL, getMinBBox)
    bbFoM   <- lapply(bbs, function(x) { x$FoM } )
    bbDiag  <- lapply(bbs, function(x) { x$diag } )
    bbFoML  <- Map(makeMOA, bbFoM,  dst=dstTarget, conversion=conversion)
    bbDiagL <- Map(makeMOA, bbDiag, dst=dstTarget, conversion=conversion)
    res$bbFoM  <- do.call("cbind", bbFoML)
    res$bbDiag <- do.call("cbind", bbDiagL)

    ## for axis limits
    bbMinPts  <- do.call("rbind", lapply(bbs, function(x) x$pts))
    axisLimsX <- c(axisLimsX, bbMinPts[ , 1])
    axisLimsY <- c(axisLimsY, bbMinPts[ , 2])

    ## radius of minimum enclosing circle
    minCircs    <- lapply(xyL, getMinCircle)
    minCircRad  <- lapply(minCircs, function(x) { x$rad } )     # radius
    minCircRadL <- Map(makeMOA, minCircRad, dst=dstTarget, conversion=conversion)
    res$minCircleRad <- do.call("cbind", minCircRadL)

    ## for axis limits
    getMinCircLims <- function(x) {
        cbind(X=c(x$ctr[1] + x$rad, x$ctr[1] - x$rad),
              Y=c(x$ctr[2] + x$rad, x$ctr[2] - x$rad))
    }
    minCircLims <- do.call("rbind", lapply(minCircs, getMinCircLims))
    axisLimsX   <- c(axisLimsX, minCircLims[ , 1])
    axisLimsY   <- c(axisLimsY, minCircLims[ , 2])

    ## Rayleigh sigma and MR
    rayParL <- Map(getRayParam, xyL, level=CIlevel, doRob=FALSE)
    sigma   <- lapply(rayParL, function(x) { x$sigma["sigma"]   })
    MR      <- lapply(rayParL, function(x) { x$MR["MR"]   })
    sigMRci <- lapply(rayParL, function(x) {
        sigMRci <- c(x$sigma[c("sigCIlo", "sigCIup")],
                     x$MR[   c("MRciLo",  "MRciUp")])
        setNames(sigMRci, c("sigma (", "sigma )", "MR (", "MR )"))
    })

    sigmaL <- Map(makeMOA, sigma, dst=dstTarget, conversion=conversion)
    MRL    <- Map(makeMOA, MR,    dst=dstTarget, conversion=conversion)
    res$sigma <- do.call("cbind", sigmaL)
    res$MR    <- do.call("cbind", MRL)
    res$sigmaMRci <- Map(makeMOA, sigMRci, dst=dstTarget, conversion=conversion)

    ## 50% circular error probable
    CEPlist <- Map(getCEP, xyL, CEPlevel=CEPlevel, dstTarget=dstTarget,
                   conversion=conversion, type=CEPtype, doRob=FALSE)
    CEPl    <- lapply(CEPlist, function(x) { x$CEP[[1]][ , CEPtype, drop=FALSE] })
    CEPmat  <- do.call("cbind", CEPl)    # as matrix
    colnames(CEPmat) <- names(xyL)
    res$CEP <- CEPmat

    #####-----------------------------------------------------------------------
    ## tests for equal spread
    ## 2 groups:   Ansari-Bradley for x- and y-coords
    ##             Kruskal-Wallis for distance to center
    ## > 2 groups: Fligner-Killeen for x- and y-coords
    ##             Wilcoxon Rank Sum (= Mann-Whitney U) for distance to center
    dstCtrGrp <- unlist(dstCtrL)           # distances to group center grouped by series
    MRGrp     <- do.call("rbind", MRL)[ , "unit"]
    MRCIGrp   <- do.call("rbind", sigMRci)
    MRDF      <- data.frame(MR=MRGrp,
                            MRciLo=MRCIGrp[ , "MR ("],
                            MRciUp=MRCIGrp[ , "MR )"],
                            series=factor(levels(DF$series)))
    names(dstCtrGrp) <- NULL

    ## create data frame with added series factor
    dstCtrDF <- data.frame(dstCtr=dstCtrGrp,
                           series=factor(rep(seq_along(levels(DF$series)), nObs),
                                         labels=levels(DF$series)))

    if(nS == 2L) {                       # compare two groups
        res$AnsariX  <- coin::ansari_test(x ~ series, alternative=ABalt,
                                          data=DF, distribution="exact")
        res$AnsariY  <- coin::ansari_test(y ~ series, alternative=ABalt,
                                          data=DF, distribution="exact")
        res$Wilcoxon <- coin::wilcox_test(dstCtr ~ series, alternative=Walt,
                                          data=dstCtrDF, distribution="exact")
    } else {                             # compare more than two groups
        res$FlignerX <- coin::fligner_test(x ~ series, data=DF,
                                           distribution=coin::approximate(B=9999))  # x
        res$FlignerY <- coin::fligner_test(y ~ series, data=DF,
                                           distribution=coin::approximate(B=9999))  # y
        res$Kruskal  <- coin::kruskal_test(dstCtr ~ series,    # dist to center
                                           data=dstCtrDF, distribution=coin::approximate(B=9999))
    }

    if(plots) {
        ## infer (x,y)-coord units from conversion
        unitXY  <- getUnits(conversion, first=FALSE)
        unitDst <- getUnits(conversion, first=TRUE)
        devNew  <- getDevice()           # platform-dependent window open

        #####-------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        syms <- c(4, 16, 2, 1, 6, 8, 3, 5, 7, 9:13, 15, 17:25)  # data symbols
        cols <- getColors(nS)            # colors

        if(nS > length(syms)) {
            stop(paste("At most", length(syms), "series possible"))
        }

        ## confidence ellipse
        confElls <- Map(getConfEll, xyL, level=CEPlevel, dst=dstTarget, conversion=conversion)

        ## adjust axis limits
        getConfEllLims <- function(x) {
            cbind(X=c(x$ctr[1] + x$size["unit", "semi-major"],
                      x$ctr[1] - x$size["unit", "semi-major"]),
                  Y=c(x$ctr[2] + x$size["unit", "semi-major"],
                      x$ctr[2] - x$size["unit", "semi-major"]))
        }

        confEllLims <- do.call("rbind", lapply(confElls, getConfEllLims))
        axisLimsX   <- c(axisLimsX, confEllLims[ , 1])
        axisLimsY   <- c(axisLimsY, confEllLims[ , 2])

        ## determine axis limits
        xLims <- range(c(DF$x, axisLimsX))
        yLims <- range(c(DF$y, axisLimsY))

        devNew()                         # open new diagram
        plot(y ~ x, data=DF, xlim=xLims, ylim=yLims, asp=1, lwd=2,
             pch=syms[unclass(DF$series)], col=cols[unclass(DF$series)],
             main=paste0("Groups with ", 100*CEPlevel, "% confidence ellipse", sep=""),
             sub=paste("distance:", dstTarget, unitDst),
             xlab=paste0("X [", unitXY, "]"), ylab=paste0("Y [", unitXY, "]"))
        abline(v=0, h=0, col="lightgray")  # add point of aim

        ## add confidence ellipses and group centers
        for(i in seq(along=xyL)) {
            drawEllipse(confElls[[i]], fg=cols[i],
                        lwd=2, pch=syms[i], cex=3)
            points(res$ctr[1, i], res$ctr[2, i], pch=syms[i], col=cols[i],
                   cex=3, lwd=2)
        }

        ## add legend
        legend(x="bottomleft", legend=names(xyL), lty=NA, pch=syms[seq_len(nS)],
               col=cols[seq_len(nS)], lwd=2, pt.cex=1.5, bg=rgb(1, 1, 1, 0.6))

        #####-------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        devNew()                         # open new diagram
        plot(y ~ x, data=DF, asp=1, xlim=xLims, ylim=yLims, lwd=2,
             pch=syms[unclass(DF$series)], col=cols[unclass(DF$series)],
             main="Groups w/ minimum bounding box & maximum spread",
             sub=paste("distance:", dstTarget, unitDst),
             xlab=paste0("X [", unitXY, "]"), ylab=paste0("Y [", unitXY, "]"))
        abline(v=0, h=0, col="lightgray")  # add point of aim
        points(res$ctr[1, ], res$ctr[2, ], col=cols[seq_len(nS)],
               pch=syms[seq_len(nS)], lwd=2, cex=3)

        ## add bounding box and maximum group spread
        for(i in seq(along=xyL)) {
            bb <- bbs[[i]]
            ## drawBox(bb, fg=cols[i])
            drawBox2(bb, fg=cols[i])
            segments(x0=xyL[[i]][maxPDidx[1, i], 1], y0=xyL[[i]][maxPDidx[1, i], 2],
                     x1=xyL[[i]][maxPDidx[2, i], 1], y1=xyL[[i]][maxPDidx[2, i], 2],
                     col=cols[i])
        }

        ## add legend
        legend(x="bottomleft", legend=names(xyL), lty=NA, pch=syms[seq_len(nS)],
               col=cols[seq_len(nS)], lwd=2, pt.cex=1.5, bg=rgb(1, 1, 1, 0.6))

        #####-------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        devNew()                         # open new diagram
        plot(y ~ x, data=DF, asp=1, xlim=xLims, ylim=yLims, lwd=2,
             pch=syms[unclass(DF$series)], col=cols[unclass(DF$series)],
             main="Groups w/ minimum enclosing circle and mean dist to center",
             sub=paste("distance:", dstTarget, unitDst),
             xlab=paste0("X [", unitXY, "]"), ylab=paste0("Y [", unitXY, "]"))
        abline(v=0, h=0, col="lightgray")  # add point of aim
        points(res$ctr[1, ], res$ctr[2, ], col=cols[seq_len(nS)], pch=syms[seq_len(nS)],
               lwd=2, cex=3)

        ## add circle with mean distance to center and minimum enclosing circle
        for(i in seq(along=xyL)) {
            drawCircle(res$ctr[ , i], radius=meanDstCtr[[i]], fg=cols[i])
            drawCircle(minCircs[[i]], fg=cols[i])
        }

        ## add legend
        legend(x="bottomleft", legend=names(xyL), lty=NA, pch=syms[seq_len(nS)],
               col=cols[seq_len(nS)], lwd=2, pt.cex=1.5, bg=rgb(1, 1, 1, 0.6))

        #####-------------------------------------------------------------------
        ## diagram: distances to center
        ## grouped boxplot + Rayleigh MR+CI
        devNew()

        op <- par(mfrow=c(1, 2))
        yLims <- c(0, max(dstCtrDF$dstCtr, MRDF$MRciUp))
        boxplot(dstCtr ~ series, data=dstCtrDF,
                main="Distance to center",
                sub=paste("distance:", dstTarget, unitDst),
                xaxt="n", col=cols, ylim=yLims,
                xlab="group", ylab=paste0("distance to center [", unitXY, "]"))
        axis(side=1, at=seq_along(levels(dstCtrDF$series)),
             labels=substring(levels(dstCtrDF$series), 1, 7), las=2)

        ## Rayleigh MR+CI
        ## raw distances to center
        stripchart(dstCtr ~ series, data=dstCtrDF, pch=20, vert=TRUE,
                   method="jitter",
                   main=paste0("Dist2ctr w/ Rayleigh MR + ", 100*CIlevel, "% CI"),
                   sub=paste("distance:", dstTarget, unitDst),
                   xaxt="n", col=adjustcolor(cols, alpha.f=0.5),
                   xlim=range(seq_along(levels(dstCtrDF$series))) + c(-0.5, 0.5),
                   ylim=yLims,
                   xlab="group", ylab=paste0("distance to center [", unitXY, "]"))
        axis(side=1, at=seq_along(levels(dstCtrDF$series)),
             labels=substring(levels(dstCtrDF$series), 1, 7), las=2)

        ## MR per group
        with(MRDF,
             points(seq_along(series), MR, pch=4, cex=2, col="black", lwd=2))
        ## MR CI per group
        with(MRDF,
             arrows(x0=seq_along(series), y0=MRciLo,
                    x1=seq_along(series), y1=MRciUp,
                    code=3, angle=90, length=0.1, col="black", lwd=2))

        par(op)
    }                                    # if(plots)

    #####-----------------------------------------------------------------------
    ## return all the collected numerical results and tests
    return(res)
}

compareGroupsPlot <-
function(DF, which=1L, xyTopLeft=TRUE, CEPlevel=0.5, CIlevel=0.95, conversion="m2cm") {
    if(!is.data.frame(DF)) { stop("DF must be a data frame") }

    which <- match.arg(as.character(which), choices=1:4)

    #####-----------------------------------------------------------------------
    ## make sure DF has the required variable names and at least two groups
    varNames <- tolower(names(DF))       # what variables are present
    needsSer <- "series"                 # required
    needsXY1 <- c("point.x", "point.y")  # coordinates must have this name
    needsXY2 <- c("x", "y")              # or this
    wantsDst <- "distance"               # useful
    wantsAIM <- c("aim.x", "aim.y")      # useful
    hasSer   <- needsSer %in% varNames   # required ones we have
    hasXY1   <- needsXY1 %in% varNames   # coordinates we have
    hasXY2   <- needsXY2 %in% varNames
    hasDst   <- wantsDst %in% varNames   # useful ones we have
    hasAIM   <- wantsAIM %in% varNames   # useful ones we have

    if(!all(hasSer)) {
        stop(c("The data frame is missing variable\n",
               paste(needsSer[!hasSer], collapse=" ")))
    }

    if(!xor(all(hasXY1), all(hasXY2))) {
        stop("Coordinates must be named either X, Y or Point.X, Point.Y")
    }

    if(!all(hasDst)) {
        warning(c("The data frame is missing variable\n",
                  paste(wantsDst[!hasDst], collapse=" "), "\n",
                  "Distance is assumed to be 100"))
        DF$distance <- 100
    }

    if(!all(hasAIM)) {
        warning(c("The data frame is missing variable(s)\n",
                  paste(wantsAIM[!hasAIM], collapse=" "), "\n",
                  "Point of Aim is assumed to be in (0,0)"))
        DF$aim.x <- 0
        DF$aim.y <- 0
    }

    #####-----------------------------------------------------------------------
    ## prepare data
    res <- vector("list", 0)               # empty list to later collect the results
    if(!is.factor(DF$series)) {            # make sure series is a factor
        DF$series <- as.factor(DF$series)
    } else {
        DF$series <- droplevels(DF$series) # remove all non-used factor levels
    }

    ## check if we have enough groups and points per group
    if(nlevels(DF$series) < 2L)            { stop("We need >= 2 groups for a comparison") }
    if(any(xtabs(~ series, data=DF) < 2L)) { stop("We need >= 2 points in each group") }

    ## prepare data: get (x,y)-coords relative to point of aim as matrix
    xy <- getXYmat(DF, xyTopLeft=xyTopLeft)
    DF <- cbind(DF, xy)
    dstTarget <- tapply(DF$distance, DF$series, mean)  # distances to target

    ## for each group extract the new (x,y)-coords as a matrix
    extractXY <- function(x) {
        DFsub <- x[ , which(!(names(x) %in% c("x", "y")))]
        getXYmat(DFsub, xyTopLeft=xyTopLeft)
    }

    xyL  <- lapply(split(DF, DF$series), extractXY)
    nS   <- length(xyL)                   # total number of series
    nObs <- vapply(xyL, nrow, integer(1)) # number of obs per series
    names(xyL) <- levels(DF$series)

    ## to determine axis limits later, collect all results in a vector
    axisLimsX <- numeric(0)
    axisLimsY <- numeric(0)

    #####-----------------------------------------------------------------------
    ## location measures
    res$ctr <- sapply(xyL, colMeans)     # group centers

    ## (mean) distances to group center
    dstCtrL    <- lapply(xyL, function(x) getDistToCtr(x))
    meanDstCtr <- lapply(dstCtrL, mean)

    ## maximum pairwise distance = maximum group spread
    maxPD     <- lapply(xyL, getMaxPairDist)   # max pairwise distance
    maxSpread <- lapply(maxPD, function(x) { x$d } )
    maxPDidx  <- sapply(maxPD, function(x) { x$idx } )

    ## bounding box figure of merit and diagonal
    bbs    <- lapply(xyL, getMinBBox)
    bbDiag <- lapply(bbs, function(x) { x$diag } )

    ## for axis limits
    bbMinPts  <- do.call("rbind", lapply(bbs, function(x) x$pts))
    axisLimsX <- c(axisLimsX, bbMinPts[ , 1])
    axisLimsY <- c(axisLimsY, bbMinPts[ , 2])

    ## radius of minimum enclosing circle
    minCircs   <- lapply(xyL, getMinCircle)
    minCircRad <- lapply(minCircs, function(x) { x$rad } )     # radius

    ## for axis limits
    getMinCircLims <- function(x) {
        cbind(X=c(x$ctr[1] + x$rad, x$ctr[1] - x$rad),
              Y=c(x$ctr[2] + x$rad, x$ctr[2] - x$rad))
    }

    minCircLims <- do.call("rbind", lapply(minCircs, getMinCircLims))
    axisLimsX   <- c(axisLimsX, minCircLims[ , 1])
    axisLimsY   <- c(axisLimsY, minCircLims[ , 2])

    ## Rayleigh sigma and MR
    rayParL <- Map(getRayParam, xyL, level=CIlevel, doRob=FALSE)
    sigma   <- lapply(rayParL, function(x) { x$sigma["sigma"]   })
    MR      <- lapply(rayParL, function(x) { x$MR["MR"]   })
    sigMRci <- lapply(rayParL, function(x) {
        sigMRci <- c(x$sigma[c("sigCIlo", "sigCIup")],
                     x$MR[   c("MRciLo",  "MRciUp")])
        setNames(sigMRci, c("sigma (", "sigma )", "MR (", "MR )"))
    })

    sigmaL <- Map(makeMOA, sigma, dst=dstTarget, conversion=conversion)
    MRL    <- Map(makeMOA, MR,    dst=dstTarget, conversion=conversion)
    res$sigma <- do.call("cbind", sigmaL)
    res$MR    <- do.call("cbind", MRL)
    res$sigmaMRci <- Map(makeMOA, sigMRci, dst=dstTarget, conversion=conversion)

    ## distance to center, Rayleigh sigma + MR
    dstCtrGrp <- unlist(dstCtrL)           # distances to group center grouped by series
    MRGrp     <- do.call("rbind", MRL)[ , "unit"]
    MRCIGrp   <- do.call("rbind", sigMRci)
    MRDF      <- data.frame(MR=MRGrp,
                            MRciLo=MRCIGrp[ , "MR ("],
                            MRciUp=MRCIGrp[ , "MR )"],
                            series=factor(levels(DF$series)))
    names(dstCtrGrp) <- NULL

    ## create data frame with added series factor
    dstCtrDF <- data.frame(dstCtr=dstCtrGrp,
                           series=factor(rep(seq_along(levels(DF$series)), nObs),
                                         labels=levels(DF$series)))

    ## confidence ellipse
    confElls <- Map(getConfEll, xyL, level=CEPlevel, dst=dstTarget, conversion=conversion)

    ## for axis limits
    getConfEllLims <- function(x) {
        cbind(X=c(x$ctr[1] + x$size["unit", "semi-major"],
                  x$ctr[1] - x$size["unit", "semi-major"]),
              Y=c(x$ctr[2] + x$size["unit", "semi-major"],
                  x$ctr[2] - x$size["unit", "semi-major"]))
    }

    confEllLims <- do.call("rbind", lapply(confElls, getConfEllLims))
    axisLimsX   <- c(axisLimsX, confEllLims[ , 1])
    axisLimsY   <- c(axisLimsY, confEllLims[ , 2])

    ## infer (x,y)-coord units from conversion
    unitXY  <- getUnits(conversion, first=FALSE)
    unitDst <- getUnits(conversion, first=TRUE)

    ## determine axis limits
    xLims <- range(c(DF$x, axisLimsX))
    yLims <- range(c(DF$y, axisLimsY))

    syms <- c(4, 16, 2, 1, 6, 8, 3, 5, 7, 9:13, 15, 17:25)  # data symbols
    cols <- getColors(nS)            # colors

    if(nS > length(syms)) {
        stop(paste("At most", length(syms), "series possible"))
    }

    if(which == 1L) {
        #####-----------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        plot(y ~ x, data=DF, xlim=xLims, ylim=yLims, asp=1, lwd=2,
             pch=syms[unclass(DF$series)], col=cols[unclass(DF$series)],
             main=paste0("Groups with ", 100*CEPlevel, "% confidence ellipse", sep=""),
             sub=paste("distance:", dstTarget, unitDst),
             xlab=paste0("X [", unitXY, "]"), ylab=paste0("Y [", unitXY, "]"))
        abline(v=0, h=0, col="lightgray")  # add point of aim

        ## add confidence ellipses and group centers
        for(i in seq(along=xyL)) {
            drawEllipse(confElls[[i]], fg=cols[i],
                        lwd=2, pch=syms[i], cex=3)
            points(res$ctr[1, i], res$ctr[2, i], pch=syms[i], col=cols[i],
                   cex=3, lwd=2)
        }

        ## add legend
        legend(x="bottomleft", legend=names(xyL), lty=NA, pch=syms[seq_len(nS)],
               col=cols[seq_len(nS)], lwd=2, pt.cex=1.5, bg=rgb(1, 1, 1, 0.6))
    }

    if(which == 2L) {
        #####-----------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        plot(y ~ x, data=DF, asp=1, xlim=xLims, ylim=yLims, lwd=2,
             pch=syms[unclass(DF$series)], col=cols[unclass(DF$series)],
             main="Groups w/ minimum bounding box & maximum spread",
             sub=paste("distance:", dstTarget, unitDst),
             xlab=paste0("X [", unitXY, "]"), ylab=paste0("Y [", unitXY, "]"))
        abline(v=0, h=0, col="lightgray")  # add point of aim
        points(res$ctr[1, ], res$ctr[2, ], col=cols[seq_len(nS)],
               pch=syms[seq_len(nS)], lwd=2, cex=3)

        ## add bounding box and maximum group spread
        for(i in seq(along=xyL)) {
            bb <- bbs[[i]]
            ## drawBox(bb, fg=cols[i])
            drawBox2(bb, fg=cols[i])
            segments(x0=xyL[[i]][maxPDidx[1, i], 1], y0=xyL[[i]][maxPDidx[1, i], 2],
                     x1=xyL[[i]][maxPDidx[2, i], 1], y1=xyL[[i]][maxPDidx[2, i], 2],
                     col=cols[i])
        }

        ## add legend
        legend(x="bottomleft", legend=names(xyL), lty=NA, pch=syms[seq_len(nS)],
               col=cols[seq_len(nS)], lwd=2, pt.cex=1.5, bg=rgb(1, 1, 1, 0.6))
    }

    if(which == 3L) {
        #####-----------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        plot(y ~ x, data=DF, asp=1, xlim=xLims, ylim=yLims, lwd=2,
             pch=syms[unclass(DF$series)], col=cols[unclass(DF$series)],
             main="Groups w/ minimum enclosing circle and mean dist to center",
             sub=paste("distance:", dstTarget, unitDst),
             xlab=paste0("X [", unitXY, "]"), ylab=paste0("Y [", unitXY, "]"))
        abline(v=0, h=0, col="lightgray")  # add point of aim
        points(res$ctr[1, ], res$ctr[2, ], col=cols[seq_len(nS)], pch=syms[seq_len(nS)],
               lwd=2, cex=3)

        ## add circle with mean distance to center and minimum enclosing circle
        for(i in seq(along=xyL)) {
            drawCircle(res$ctr[ , i], radius=meanDstCtr[[i]], fg=cols[i])
            drawCircle(minCircs[[i]], fg=cols[i])
        }

        ## add legend
        legend(x="bottomleft", legend=names(xyL), lty=NA, pch=syms[seq_len(nS)],
               col=cols[seq_len(nS)], lwd=2, pt.cex=1.5, bg=rgb(1, 1, 1, 0.6))
    }                                    # if(plots)

    if(which == 4L) {
        #####-------------------------------------------------------------------
        ## diagram: distances to center
        ## grouped boxplot + Rayleigh MR+CI
        op <- par(mfrow=c(1, 2))
        yLims <- c(0, max(dstCtrDF$dstCtr, MRDF$MRciUp))
        boxplot(dstCtr ~ series, data=dstCtrDF,
                main="Distance to center",
                sub=paste("distance:", dstTarget, unitDst),
                xaxt="n", col=cols, ylim=yLims,
                xlab="group", ylab=paste0("distance to center [", unitXY, "]"))
        axis(side=1, at=seq_along(levels(dstCtrDF$series)),
             labels=substring(levels(dstCtrDF$series), 1, 7), las=2)

        ## Rayleigh MR+CI
        stripchart(dstCtr ~ series, data=dstCtrDF, pch=20, vert=TRUE,
                   method="jitter",
                   main=paste0("Dist2ctr w/ Rayleigh MR + ", 100*CIlevel, "% CI"),
                   sub=paste("distance:", dstTarget, unitDst),
                   xaxt="n", col=adjustcolor(cols, alpha.f=0.5),
                   xlim=range(seq_along(levels(dstCtrDF$series))) + c(-0.5, 0.5),
                   ylim=yLims,
                   xlab="group", ylab=paste0("distance to center [", unitXY, "]"))
        axis(side=1, at=seq_along(levels(dstCtrDF$series)),
             labels=substring(levels(dstCtrDF$series), 1, 7), las=2)

        with(MRDF,
             points(seq_along(series), MR, pch=4, cex=2, col="black", lwd=2))
        with(MRDF,
             arrows(x0=seq_along(series), y0=MRciLo,
                    x1=seq_along(series), y1=MRciUp,
                    code=3, angle=90, length=0.1, col="black", lwd=2))
        par(op)
    }

    #####-----------------------------------------------------------------------
    ## return all the collected numerical results and tests
    return(res)
}
