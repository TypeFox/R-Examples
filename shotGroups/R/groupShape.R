groupShape <-
function(xy, plots=TRUE, bandW=0.5, outlier=c("mcd", "pca"),
         dstTarget=100, conversion="m2cm", ...) {
    UseMethod("groupShape")
}

groupShape.data.frame <-
function(xy, plots=TRUE, bandW=0.5, outlier=c("mcd", "pca"),
         dstTarget=100, conversion="m2cm", ...) {
    xy <- getXYmat(xy)
    NextMethod("groupShape")
}

groupShape.default <-
function(xy, plots=TRUE, bandW=0.5, outlier=c("mcd", "pca"),
         dstTarget=100, conversion="m2cm", ...) {
    if(!is.matrix(xy))  { stop("xy must be a matrix") }
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(ncol(xy) != 2L)  { stop("xy must have two columns") }

    outlier <- match.arg(outlier)

    #####-----------------------------------------------------------------------
    ## prepare data
    X      <- xy[ , 1]                   # x-coords
    Y      <- xy[ , 2]                   # y-coords
    Npts   <- nrow(xy)                   # number of points
    res    <- vector("list", 0)          # empty list to later collect the results
    devNew <- getDevice()                # platform-dependent window open

    haveMVoutlier <- requireNamespace("mvoutlier", quietly=TRUE)
    haveEnergy    <- requireNamespace("energy",    quietly=TRUE)

    haveRob <- if(Npts < 4L) {
        warning(c("We need >= 4 points for robust estimations,\n",
                  "outlier analysis, and chi^2 plot for multivariate normality"))
        FALSE                           # can we do robust estimation?
    }  else {
        TRUE
    }                                   # if(Npts < 4L)

    #####-----------------------------------------------------------------------
    ## correlation matrix of (x,y)-coords
    res$corXY <- cor(xy)

    if(haveRob) {
        rob <- robustbase::covMcd(xy, cor=TRUE)
        res$corXYrob <- rob$cor         # robust estimate correlation-matrix

        #####-------------------------------------------------------------------
        ## outlier-analysis for joint distribution of (x,y)-coords
        if(plots && haveMVoutlier) {
            devNew()                         # open new diagram
            op <- par(no.readonly=TRUE)      # save device parameters

            ## outlier-analysis-plot         # this can fail due to memory constraints
            outXY <- tryCatch(mvoutlier::aq.plot(xy), error=function(e) {
                warning(c("mvoutlier::aq.plot() failed:\n", e$message))
                return(list(outliers=NA)) })
            par(op)                          # reset device parameters
        } else if(haveMVoutlier) {
            pdf(file=NULL)                   # redirect diagram to null device
            ## outlier-analysis-plot         # this can fail due to memory constraints
            outXY <- tryCatch(mvoutlier::aq.plot(xy), error=function(e) {
                warning(c("mvoutlier::aq.plot() failed:\n", e$message))
                return(list(outliers=NA)) })
            dev.off()
        } else {
            warning(c("Package mvoutlier is not installed"))
            outXY <- list(outliers=NA)
        }

        ##  identified outliers
        res$Outliers <- if((outlier == "mcd") && !is.na(outXY$outliers)) {
            which(outXY$outliers)
        } else if((outlier == "pca") && haveMVoutlier) {
            which(mvoutlier::pcout(xy, makeplot=FALSE, ...)$wfinal01 == 0)
        } else { NULL }
    } else {
        res$corXYrob <- NULL
        res$Outliers <- NULL
    }                                    # if(haveRob)

    #####-----------------------------------------------------------------------
    ## normality tests
    ## Shapiro-Wilk-Tests for normality of (x,y)-coords separately
    if(Npts >= 3L) {
        if(Npts <= 5000) {               # Shapiro-Wilk only for <= 5000 points
            res$ShapiroX <- shapiro.test(X)  # normality x-coords
            res$ShapiroY <- shapiro.test(Y)  # normality y-coords
        } else {
            warning(c("Shapiro-Wilk normality test only supports <= 5000 points\n",
                      "will use plug-in Kolmogorov-Smirnov-Test instead"))
            ## Gauss bias correction factor for estimated sd
            Ncorr   <- 2*(Npts-1)            # effective N
            corrFac <- 1 / (sqrt(2/(Ncorr-1)) * exp(lgamma(Ncorr/2) - lgamma((Ncorr-1)/2)))

            corrG   <- ifelse(is.finite(corrFac), corrFac, 1)
            sigHatX <- corrG * sd(X)
            sigHatY <- corrG * sd(Y)

            ## Kolmogorov-Smirnov-Test with estimated mean and sigma
            ## this plug-in test is not strictly valid (Lilliefors)
            res$ksX <- ks.test(X, "pnorm", mean=mean(X), sd=sigHatX)
            res$ksY <- ks.test(Y, "pnorm", mean=mean(Y), sd=sigHatY)
        }

        ## Energy-Test for multivariate normality of joint (x,y)-distribution
        res$multNorm <- if(haveEnergy) {
            energy::mvnorm.etest(xy, R=1499)
        } else {
            warning("Package energy for multivariate normality test not installed")
            NULL
        }
   } else {
        res$ShapiroX <- NULL
        res$ShapiroY <- NULL
        res$multNorm <- NULL
        warning("We need >= 3 points for normality tests")
    }                                    # if(Npts >= 3L)

    if(plots) {
        ## infer (x,y)-coord units from conversion
        unitXY  <- getUnits(conversion, first=FALSE)
        unitDst <- getUnits(conversion, first=TRUE)

        ## to determine axis limits later, collect all results in a vector
        axisLimsX <- numeric(0)
        axisLimsY <- numeric(0)

        #####-------------------------------------------------------------------
        ## diagram: separate Q-Q-plots for eyeballing normality in x- and y-coords
        devNew()                         # open new diagram
        qqnorm(X, pch=20, main="Q-Q-plot x-coordinates for eyeballing normality",
               sub=paste("distance:", dstTarget, unitDst),
               xlab="Quantiles from standard normal distribution",
               ylab=paste0("Observed quantiles [", unitXY, "]"))
        qqline(X, col="red", lwd=2)      # reference line

        devNew()                         # open new diagram
        qqnorm(Y, pch=20, main="Q-Q-plot y-coordinates for eyeballing normality",
               sub=paste("distance:", dstTarget, unitDst),
               xlab="Quantiles from standard normal distribution",
               ylab=paste0("Observed quantiles [", unitXY, "]"))
        qqline(Y, col="red", lwd=2)      # reference line

        #####-------------------------------------------------------------------
        ## diagram: histograms for x- and y-coords
        ## x-coords
        ## choose y-axis limits
        maxNormX <- getMaxNorm(X, 2)[2]
        hX       <- hist(X, breaks="FD", plot=FALSE)
        yLims    <- c(0, max(c(hX$density, maxNormX)))

        devNew()                         # open new diagram
        plot(hX, ylim=yLims, freq=FALSE,
             main="Histogram x-coordinates w/ kernel density estimate",
             sub=paste("distance:", dstTarget, unitDst),
             xlab=paste0("X [", unitXY, "]"))
        rug(jitter(X))                   # show single values

        ## add fitted normal curve and kernel density estimate
        dnormX <- function(x) { dnorm(x, mean(X), sd(X)) }
        curve(dnormX, lwd=2, col="blue", add=TRUE)
        lines(density(X), lwd=2, col="red")  # kernel density estimate

        ## add legend
        legend(x="topleft", legend=c("normal distribution", "kernel density estimate"),
               col=c("blue", "red"), lty=c(1, 1), lwd=c(2, 2), bg=rgb(1, 1, 1, 0.7))

        ## histogram y-coords
        ## choose y-axis limits
        maxNormY <- getMaxNorm(Y, 2)[2]
        hY       <- hist(Y, breaks="FD", plot=FALSE)
        yLims    <- c(0, max(c(hY$density, maxNormY)))

        devNew()                         # open new diagram
        plot(hY, ylim=yLims, freq=FALSE,
             main="Histogram y-coordinates w/ kernel density estimate",
             sub=paste("distance:", dstTarget, unitDst),
             xlab=paste0("Y [", unitXY, "]"))
        rug(jitter(Y))                   # show single values

        ## add fitted normal curve and kernel density estimate
        dnormY <- function(x) { dnorm(x, mean(Y), sd(Y)) }
        curve(dnormY, lwd=2, col="blue", add=TRUE)
        lines(density(Y), lwd=2, col="red")  # kernel density estimate

        ## add legend
        legend(x="topleft",
               legend=c("normal distribution", "kernel density estimate"),
               col=c("blue", "red"), lty=c(1, 1), lwd=c(2, 2), bg=rgb(1, 1, 1, 0.7))

        ## chi-square qq-plot for eyeballing multivariate normality
        ## quantiles of squared robust Mahalanobis distance against quantiles
        ## from chi^2 distribution with 2 df
        if(haveRob) {
            ctrRob   <- rob$center       # robust estimate: group center,
            covXYrob <- rob$cov          # group covariance matrix

            ## for axis limits
            ellSize   <- 2*sqrt(eigen(covXYrob)$values)
            axisLimsX <- c(axisLimsX, ctrRob[1] + ellSize[1],
                                      ctrRob[1] - ellSize[1])
            axisLimsY <- c(axisLimsY, ctrRob[2] + ellSize[1],
                                      ctrRob[2] - ellSize[1])

            ## squared robust Mahalanobis-distance
            mDstSq <- mahalanobis(xy, center=ctrRob, cov=covXYrob)
            devNew()                     # open new diagram
            plot(qchisq(ppoints(mDstSq), df=ncol(xy)), sort(mDstSq),
                 xlab="Quantiles chi^2 distribution",
                 ylab="Quantiles (robust Mahalanobis distances)^2", pch=20,
                 main="Chi^2 Q-Q-plot for eyeballing multivariate normality")
            abline(a=0, b=1, col="red", lwd=2)  # add a reference line
        }                                # if(haveRob)

        #####-------------------------------------------------------------------
        ## diagram: 2D-kernel density estimate for joint (x,y)-distribution
        ## determine axis limits
        xLims <- range(c(X, axisLimsX))
        yLims <- range(c(Y, axisLimsY))

        devNew()                         # open new diagram
        smoothScatter(X, Y, asp=1, bandwidth=bandW, xlim=xLims, ylim=yLims,
                      main="2D-kernel density estimate and error ellipses",
                      sub=paste("distance:", dstTarget, unitDst),
                      xlab=paste0("X [", unitXY, "]"), ylab=paste0("Y [", unitXY, "]"))
        abline(h=0, v=0, lwd=2)          # add point of aim

        ## add group center and robust estimate for group center
        if(haveRob) {
            points(ctrRob[1], ctrRob[2], col=rgb(0.5, 0.5, 0.5, 0.4),
                   pch=4, lwd=2, cex=1.5)

            ## add robust error ellipses with radius = 1 and = 2
            drawEllipse(ctrRob, covXYrob, radius=1, fg=rgb(0.5, 0.5, 0.5, 0.4),
                        pch=4, lwd=2)
            drawEllipse(ctrRob, covXYrob, radius=2, fg=rgb(0.5, 0.5, 0.5, 0.4),
                        pch=4, lwd=2)
        }                                # if(haveRob)

        ## add legend
        legend(x="bottomleft", legend=c("robust center", "robust error ellipses"),
               col=c("darkgray", "darkgray"), pch=c(4, NA), lty=c(NA, 1),
               lwd=c(2, 2), bg=rgb(1, 1, 1, 0.7))
    }                                    # if(plots)

    #####-----------------------------------------------------------------------
    ## return all the collected numerical results and tests
    return(res)
}

groupShapePlot <-
function(xy, which=1L, bandW=0.5, outlier=c("mcd", "pca"),
         dstTarget=100, conversion="m2cm", ...) {

    if(!is.data.frame(xy)) { stop("xy must be a data.frame") }
    xy <- getXYmat(xy)
    if(!is.numeric(xy))    { stop("xy must be numeric") }
    if(ncol(xy) != 2L)     { stop("xy must have two columns") }

    which   <- match.arg(as.character(which), choices=1:7)
    outlier <- match.arg(outlier)

    #####-----------------------------------------------------------------------
    ## prepare data
    X    <- xy[ , 1]                     # x-coords
    Y    <- xy[ , 2]                     # y-coords
    Npts <- nrow(xy)                     # number of points
    res  <- vector("list", 0)            # empty list to later collect the results

    haveMVoutlier <- requireNamespace("mvoutlier", quietly=TRUE)

    haveRob <- TRUE                      # can we do robust estimation?
    if(Npts < 4L) {
        warning(c("We need >= 4 points for robust estimations,\n",
                  "outlier analysis, and chi^2 plot for multivariate normality"))
        haveRob <- FALSE
    } else {
        rob      <- robustbase::covMcd(xy, cor=TRUE)
        ctrRob   <- rob$center       # robust estimate: group center,
        covXYrob <- rob$cov          # group covariance matrix
    }                                # if(Npts < 4L)

    if((which == 1L) && haveRob && haveMVoutlier) {
        #####-------------------------------------------------------------------
        ## outlier-analysis for joint distribution of (x,y)-coords
        ## this can fail due to memory constraints
        op <- par(no.readonly=TRUE)      # save device parameters
        outXY <- tryCatch(mvoutlier::aq.plot(xy), error=function(e) {
            warning(c("mvoutlier::aq.plot() failed:\n", e$message))
            return(list(outliers=NA)) })
        par(op)                          # reset device parameters
    }                                    # if((which == 1) && haveRob)

    ## infer (x,y)-coord units from conversion
    unitXY  <- getUnits(conversion, first=FALSE)
    unitDst <- getUnits(conversion, first=TRUE)

    ## to determine axis limits later, collect all results in a vector
    axisLimsX <- numeric(0)
    axisLimsY <- numeric(0)

    if(which == 2L) {
        #####-------------------------------------------------------------------
        ## diagram: separate Q-Q-plots for eyeballing normality in x- and y-coords
        qqnorm(X, pch=20, main="Q-Q-plot x-coordinates for eyeballing normality",
               sub=paste("distance:", dstTarget, unitDst),
               xlab="Quantiles from standard normal distribution",
               ylab=paste0("Observed quantiles [", unitXY, "]"))
        qqline(X, col="red", lwd=2)      # reference line
    }

    if(which == 3L) {
        qqnorm(Y, pch=20, main="Q-Q-plot y-coordinates for eyeballing normality",
               sub=paste("distance:", dstTarget, unitDst),
               xlab="Quantiles from standard normal distribution",
               ylab=paste0("Observed quantiles [", unitXY, "]"))
        qqline(Y, col="red", lwd=2)      # reference line
    }

    if(which == 4L) {
        #####-------------------------------------------------------------------
        ## diagram: histograms for x- and y-coords
        ## x-coords
        ## choose y-axis limits
        maxNormX <- getMaxNorm(X, 2)[2]
        hX       <- hist(X, breaks="FD", plot=FALSE)
        yLims    <- c(0, max(c(hX$density, maxNormX)))

        plot(hX, ylim=yLims, freq=FALSE,
             main="Histogram x-coordinates w/ kernel density estimate",
             sub=paste("distance:", dstTarget, unitDst),
             xlab=paste0("X [", unitXY, "]"))
        rug(jitter(X))                   # show single values

        ## add fitted normal curve and kernel density estimate
        dnormX <- function(x) { dnorm(x, mean(X), sd(X)) }
        curve(dnormX, lwd=2, col="blue", add=TRUE)
        lines(density(X), lwd=2, col="red")  # kernel density estimate

        ## add legend
        legend(x="topleft", legend=c("normal distribution", "kernel density estimate"),
               col=c("blue", "red"), lty=c(1, 1), lwd=c(2, 2), bg=rgb(1, 1, 1, 0.7))
    }

    if(which == 5L) {
        ## histogram y-coords
        ## choose y-axis limits
        maxNormY <- getMaxNorm(Y, 2)[2]
        hY       <- hist(Y, breaks="FD", plot=FALSE)
        yLims    <- c(0, max(c(hY$density, maxNormY)))

        plot(hY, ylim=yLims, freq=FALSE,
             main="Histogram y-coordinates w/ kernel density estimate",
             sub=paste("distance:", dstTarget, unitDst),
             xlab=paste0("Y [", unitXY, "]"))
        rug(jitter(Y))                   # show single values

        ## add fitted normal curve and kernel density estimate
        dnormY <- function(x) { dnorm(x, mean(Y), sd(Y)) }
        curve(dnormY, lwd=2, col="blue", add=TRUE)
        lines(density(Y), lwd=2, col="red")  # kernel density estimate

        ## add legend
        legend(x="topleft",
               legend=c("normal distribution", "kernel density estimate"),
               col=c("blue", "red"), lty=c(1, 1), lwd=c(2, 2), bg=rgb(1, 1, 1, 0.7))
    }

    if(which == 6L) {
        ## chi-square qq-plot for eyeballing multivariate normality
        ## quantiles of squared robust Mahalanobis distance against quantiles
        ## from chi^2 distribution with 2 df
        if(haveRob) {
            ## for axis limits
            ellSize   <- 2*sqrt(eigen(covXYrob)$values)
            axisLimsX <- c(axisLimsX, ctrRob[1] + ellSize[1],
                                      ctrRob[1] - ellSize[1])
            axisLimsY <- c(axisLimsY, ctrRob[2] + ellSize[1],
                                      ctrRob[2] - ellSize[1])

            ## squared robust Mahalanobis-distance
            mDstSq <- mahalanobis(xy, center=ctrRob, cov=covXYrob)
            plot(qchisq(ppoints(mDstSq), df=ncol(xy)), sort(mDstSq),
                 xlab="Quantiles chi^2 distribution",
                 ylab="Quantiles (robust Mahalanobis distances)^2", pch=20,
                 main="Chi^2 Q-Q-plot for eyeballing multivariate normality")
            abline(a=0, b=1, col="red", lwd=2)  # add a reference line
        }                                # if(haveRob)
    }

    if(which == 7L) {
        #####-------------------------------------------------------------------
        ## diagram: 2D-kernel density estimate for joint (x,y)-distribution
        ## determine axis limits
        xLims <- range(c(X, axisLimsX))
        yLims <- range(c(Y, axisLimsY))

        smoothScatter(X, Y, asp=1, bandwidth=bandW, xlim=xLims, ylim=yLims,
                      main="2D-kernel density estimate and error ellipses",
                      sub=paste("distance:", dstTarget, unitDst),
                      xlab=paste0("X [", unitXY, "]"), ylab=paste0("Y [", unitXY, "]"))
        abline(h=0, v=0, lwd=2)          # add point of aim

        ## add group center and robust estimate for group center
        if(haveRob) {
            points(ctrRob[1], ctrRob[2], col=rgb(0.5, 0.5, 0.5, 0.4),
                   pch=4, lwd=2, cex=1.5)

            ## add robust error ellipses with radius = 1 and = 2
            drawEllipse(ctrRob, covXYrob, radius=1, fg=rgb(0.5, 0.5, 0.5, 0.4),
                        pch=4, lwd=2)
            drawEllipse(ctrRob, covXYrob, radius=2, fg=rgb(0.5, 0.5, 0.5, 0.4),
                        pch=4, lwd=2)
        }                                # if(haveRob)

        ## add legend
        legend(x="bottomleft", legend=c("robust center", "robust error ellipses"),
               col=c("darkgray", "darkgray"), pch=c(4, NA), lty=c(NA, 1),
               lwd=c(2, 2), bg=rgb(1, 1, 1, 0.7))
    }                                    # if(plots)

    #####-----------------------------------------------------------------------
    ## return all the collected numerical results and tests
    return(invisible(NULL))
}
