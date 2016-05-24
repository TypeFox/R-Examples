groupLocation <-
function(xy, level=0.95, plots=TRUE, bootCI=c("basic", "bca"),
         dstTarget=100, conversion="m2cm") {
    UseMethod("groupLocation")
}

groupLocation.data.frame <-
function(xy, level=0.95, plots=TRUE, bootCI=c("basic", "bca"),
         dstTarget=100, conversion="m2cm") {
    xy <- getXYmat(xy)
    NextMethod("groupLocation")
}

groupLocation.default <-
function(xy, level=0.95, plots=TRUE, bootCI=c("basic", "bca"),
         dstTarget=100, conversion="m2cm") {
    if(!is.matrix(xy))       { stop("xy must be a matrix") }
    if(!is.numeric(xy))      { stop("xy must be numeric") }
    if(ncol(xy) != 2L)       { stop("xy must have two columns") }
    if(!is.numeric(level))   { stop("level must be numeric") }
    if(level <= 0)           { stop("level must be > 0") }

    bootCI <- match.arg(bootCI, choices=c("none", "norm", "basic", "perc", "bca"), several.ok=TRUE)

    ## check if CI level is given in percent
    if(level >= 1) {
        while(level >= 1) { level <- level / 100 }
        warning(c("level must be in (0,1) and was set to ", level))
    }

    #####-----------------------------------------------------------------------
    ## prepare data
    X    <- xy[ , 1]                     # x-coords
    Y    <- xy[ , 2]                     # y-coords
    Npts <- nrow(xy)                     # number of observations
    res  <- vector("list", 0)            # empty list to later collect the results

    haveRob <- if(Npts < 4L) {           # can we do robust estimation?
        warning("We need >= 4 points for robust estimations")
        haveRob <- FALSE
    } else {
        TRUE
    }                                    # if(nrow(xy) < 4L)

    #####-----------------------------------------------------------------------
    ## location measures
    res$ctr <- colMeans(xy)              # center of joint (x,y)-distribution

    ## robust estimation of center
    res$ctrRob <- if(haveRob) {
        robustbase::covMcd(xy)$center
    } else {
        NULL
    }                                    # if(haveRob)

    distPOA     <- sqrt(sum(res$ctr^2))  # distance to point of aim
    res$distPOA <- makeMOA(distPOA, dst=dstTarget, conversion=conversion)

    res$distPOArob <- if(haveRob) {      # rob distance to point of aim
        distPOArob <- sqrt(sum(res$ctrRob^2))
        makeMOA(distPOArob, dst=dstTarget, conversion=conversion)
    } else {
        NULL
    }                                    # if(haveRob)

    ## Hotelling's T^2 test for equality of (x,y)-center with point of aim (0,0)
    res$Hotelling <- if(Npts > 2L) {
        anova(lm(cbind(X, Y) ~ 1), test="Hotelling-Lawley")
    } else {
        warning("We need >= 3 points for Hotelling's T^2 test")
        NULL
    }                                    # if(Npts > 2L)

    #####-----------------------------------------------------------------------
    ## confidence intervals for x- and y-coords
    ## parametric: t-CI
    alpha  <- 1-level                    # alpha-level
    tCrit  <- qt(c(alpha/2, 1-alpha/2), Npts-1)  # critical t-values left and right
    Mx     <- mean(X)                    # mean x-coords
    My     <- mean(Y)                    # mean y-coords
    sMx    <- sd(X) / sqrt(Npts)         # standard error of the mean x
    sMy    <- sd(Y) / sqrt(Npts)         # standard error of the mean y
    ctrXci <- rbind(t=rev(Mx-tCrit*sMx)) # t-CI x-coords
    ctrYci <- rbind(t=rev(My-tCrit*sMy)) # t-CI y-coords

    ## non-parametric: bootstrap-CIs for center (basic and BCa)
    if(!("none" %in% bootCI)) {          # do bootstrap CIs
        NrplMin <- 1499L                 # minimum number of replications
        Nrpl <- if("bca" %in% bootCI) {  # number of replications
            max(NrplMin, Npts+1)         # BCa needs at least number of points
        } else {
            NrplMin
        }

        ## group center for one replication
        getCtr  <- function(x, idx) { colMeans(x[idx, ]) }
        bs      <- boot::boot(xy, statistic=getCtr, R=Nrpl)  # bootstrap centers
        xCIboot <- boot::boot.ci(bs, conf=level, type=bootCI, index=1) # x
        yCIboot <- boot::boot.ci(bs, conf=level, type=bootCI, index=2) # y

        ## CI type names in output structure of boot.ci()
        CInames <- c(basic="basic", norm="normal", perc="percent", bca="bca")
        CItype  <- CInames[bootCI]
        xCImat  <- sapply(CItype, function(x) {
            len <- length(xCIboot[[x]])
            xCIboot[[x]][(len-1):len] })
        yCImat  <- sapply(CItype, function(x) {
            len <- length(yCIboot[[x]])
            yCIboot[[x]][(len-1):len] })

        ## add bootstrap CIs to parametric CI
        ctrXci <- rbind(ctrXci, t(xCImat))
        ctrYci <- rbind(ctrYci, t(yCImat))
    }

    res$ctrXci <- ctrXci
    res$ctrYci <- ctrYci
    colnames(res$ctrXci) <- c("x (", "x )")
    colnames(res$ctrYci) <- c("y (", "y )")

    if(plots) {
        ## infer (x,y)-coord units from conversion
        unitXY  <- getUnits(conversion, first=FALSE)
        unitDst <- getUnits(conversion, first=TRUE)
        devNew  <- getDevice()           # platform-dependent window open

        #####-------------------------------------------------------------------
        ## diagram: 2D-scatter plot for the (x,y)-distribution
        devNew()                         # open new diagram
        plot(Y ~ X, asp=1, main="Group (x,y)-coordinates", pch=16,
             sub=paste("distance:", dstTarget, unitDst),
             xlab=paste0("X [", unitXY, "]"), ylab=paste0("Y [", unitXY, "]"))
        abline(v=0, h=0, col="gray")     # add point of aim

        ## add (robust) group center
        points(res$ctr[1], res$ctr[2], col="blue", pch=4, lwd=2, cex=1.5)

        if(haveRob) {
            points(res$ctrRob[1], res$ctrRob[2], col="red",
                   pch=4, lwd=2, cex=1.5)
        }                                # if(haveRob)

        ## add legend
        legend(x="bottomleft", legend=c("group center", "robust group center"),
               col=c("blue", "red"), pch=4, lty=NA, lwd=2, bg=rgb(1, 1, 1, 0.7))
    }                                    # if(plots)

    #####-----------------------------------------------------------------------
    ## return all the collected numerical results and tests
    return(res)
}
