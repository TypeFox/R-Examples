plot.weibreg <- function(x,
                         fn = c("haz", "cum", "den", "sur"),
                         main = NULL,
                         xlim = NULL,
                         ylim = NULL,
                         xlab = NULL,
                         ylab = NULL,
                         new.data = x$means,
                         ...){
    if (!inherits(x, "weibreg")) stop("Works only with 'weibreg' objects.")
    ##if (x$pfixed) stop("True exponential hazards are not plotted")
    if (!(all(fn %in% c("haz", "cum", "den", "sur"))))
        stop(paste(fn, "is an illegal value of 'fn'"))

    if (length(fn) == 4){
        oldpar <- par(mfrow = c(2, 2))
        on.exit(par(oldpar))
    }
    ncov <- length(x$means)
    ns <- x$n.strata
    if (x$pfixed){
        p <- rep(x$shape, ns)
        lambda <- exp(x$coefficients[ncov + (1:ns)])
    }else{
        p <- exp(x$coefficients[ncov + (1:ns) * 2])
        lambda <- exp(x$coefficients[ncov + (1:ns) * 2 - 1])
    }
    if (ncov){
        uppe <- exp(-sum(new.data[1:ncov] * x$coefficients[1:ncov]) / p)
        lambda <- lambda * uppe
    }
    if (is.null(xlim))
        xlim <- c(min(x$y[, 1]), max(x$y[, 2]))

    npts <- 999
    xx <- seq(xlim[1], xlim[2], length = npts)
    ##if (xx[1] <= 0) xx[1] <- 0.001


    ## hazard
    if ("haz" %in% fn){
        haz <- matrix(ncol = npts, nrow = ns)
        for (i in 1:ns){
            haz[i, ] <- hweibull(xx, scale = lambda[i], shape = p[i])
        }

        if (is.null(ylim)) ylim <- c(0, max(haz))
        if (min(p) < 1) ylim[2] <- min(ylim[2], max(haz[, -1]))

        if (is.null(xlab)) xlab <- "Duration"
        if (is.null(ylab)) ylab <- "Hazard"
        if (is.null(main)) main <- "Weibull hazard function"
        plot(xx, haz[1, ], type = "l", xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, main = main, ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, haz[i, ], type = "l", lty = i)
            }
        }
        abline(h = 0)
        abline(v = 0)
    }
    ## Cumulative hazard
    if ("cum" %in% fn){

        Haz <- matrix(ncol = npts, nrow = ns)

    ##if (is.null(ylim))
        for (i in 1:ns){
            Haz[i, ] <- Hweibull(xx, scale = lambda[i], shape = p[i])
        }
        ylim <- c(0, max(Haz))
        ##if (is.null(xlab))
        xlab <- "Duration"
        ##if (is.null(ylab))
        ylab <- "Cumulative Hazard"
        ##if (is.null(main))
        main <- "Weibull cumulative hazard function"
        plot(xx, Haz[1, ], type = "l", xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, main = main, ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, Haz[i, ], type = "l", lty = i)
            }
        }
        abline(h = 0)
        abline(v = 0)
    }
    ## density
    if ("den" %in% fn){

        den <- matrix(ncol = npts, nrow = ns)
        for (i in 1:ns){
            den[i, ] <- dweibull(xx, scale = lambda[i], shape = p[i])
        }

        ##if (is.null(ylim))
        ylim <- c(0, max(den))

        if (min(p) < 1) ylim[2] <- min(max(den[, -1]))

        ##if (is.null(xlab))
        xlab <- "Duration"
        ##if (is.null(ylab))
        ylab <- "Density"
        ##if (is.null(main))
        main <- "Weibull density function"
        plot(xx, den[1, ], type = "l", xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, main = main, ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, den[i, ], type = "l", lty = i)
            }
        }
        abline(h = 0)
        abline(v = 0)
    }
    ## Survivor function
    if ("sur" %in% fn){


        sur <- matrix(ncol = npts, nrow = ns)
        for (i in 1:ns){
            sur[i, ] <- pweibull(xx, scale = lambda[i], shape = p[i],
                                 lower.tail = FALSE)
        }

        ##if (is.null(ylim))
        ylim <- c(0, 1)

        ##if (is.null(xlab))
        xlab <- "Duration"
        ##if (is.null(ylab))
        ylab <- "Survival"
        ##if (is.null(main))
        main <- "Weibull survivor function"
        plot(xx, sur[1, ], type = "l", xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, main = main, ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, sur[i, ], type = "l", lty = i)
            }
        }
        abline(h = 0)
        abline(v = 0)
    }
    ##par(oldpar)
}
