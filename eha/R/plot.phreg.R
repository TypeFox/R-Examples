plot.phreg <- function(x,
                       fn = c("haz", "cum", "den", "sur"),
                       main = NULL,
                       xlim = NULL,
                       ylim = NULL,
                       xlab = "Duration",
                       ylab = "",
                       col,   ## New 6 Feb 2013
                       lty,   ## New 6 Feb 2013
                       printLegend = TRUE,
                       ##legend = printLegend,
                       new.data = NULL,
                         ...){

    if (!inherits(x, "phreg")) stop("Works only with 'phreg' objects.")
    if (!is.null(new.data)) warning("argument 'newdata' is not used any more")
    if (missing(col)) col <- rep(1, x$n.strata) ## New 2013-12-05
    if (missing(lty)) lty <- 1:x$n.strata # No. of strata

    if (length(col) < x$n.strata) col <- rep(col, x$n.strata)
    if (length(lty) < x$n.strata) lty <- rep(lty, x$n.strata)

    if (!(all(fn %in% c("haz", "cum", "den", "sur"))))
        stop(paste(fn, "is an illegal value of 'fn'"))


    if (length(fn) >= 3){
        oldpar <- par(mfrow = c(2, 2))
        on.exit(par(oldpar))
    }else if (length(fn) == 2){
        oldpar <- par(mfrow = c(2, 1))
        on.exit(par(oldpar))
    }
    ncov <- length(x$means)
    ns <- x$n.strata
    if (x$pfixed){
        shape <- rep(x$shape, ns)
        scale <- exp(x$coefficients[ncov + (1:ns)])
    }else if (x$dist != "pch"){
        shape <- exp(x$coefficients[ncov + (1:ns) * 2])
        scale <- exp(x$coefficients[ncov + (1:ns) * 2 - 1])
    }

    if (ncov && x$center){# New in 2.4-0:
        score <- exp(sum(x$means * x$coefficients[1:ncov]))
    }else{
        score <- 1
    }

    ##if (ncov){ # THIS IS for aftplot!!
    ##    uppe <- exp(-sum(new.data[1:ncov] * x$coefficients[1:ncov]) / p)
    ##    lambda <- lambda * uppe
    ##}
    if (is.null(xlim)){
        xlim <- c(min(x$y[, 1]), max(x$y[, 2]))
        if (xlim[1] <= 0) xlim[1] <- 0.001 * xlim[2] ## Avoid Inf at 0! (hack!?)
    }
    npts <- 4999
    xx <- seq(xlim[1], xlim[2], length = npts)
    haz <- matrix(0, ncol = npts, nrow = ns)
    sur <- haz
    Haz <- haz
    ## hazard
    if (x$dist == "weibull"){
        for (i in 1:ns){
            scal <- scale[i] / (score)^(1/shape[i])
            haz[i, ] <- hweibull(xx, shape = shape[i],
                                 scale = scal)
            sur[i, ] <- pweibull(xx, shape = shape[i],
                                 scale = scal, lower.tail = FALSE)
            Haz[i, ] <- Hweibull(xx, shape = shape[i],
                                 scale = scal)
        }
        dist <- "Weibull"
    }else if (x$dist == "pch"){
        for (i in 1:ns){
            haz[i, ] <- hpch(xx, x$cuts, score * x$hazards[i, ])
            sur[i, ] <- ppch(xx, x$cuts, score * x$hazards[i, ],
                             lower.tail = FALSE)
            Haz[i, ] <- Hpch(xx, x$cuts, score * x$hazards[i, ])
        }
        dist <- "Pcwise const"
    }else if (x$dist == "loglogistic"){
        for (i in 1:ns){
            haz[i, ] <- hllogis(xx, shape = shape[i],
                                 scale = scale[i]) * score
            sur[i, ] <- pllogis(xx, shape = shape[i],
                                 scale = scale[i], lower.tail = FALSE)^score
            Haz[i, ] <- Hllogis(xx, shape = shape[i],
                                 scale = scale[i]) * score
        }
        dist <- "Loglogistic"

    }else if (x$dist == "lognormal"){
        sdlog <- 1 / shape
        meanlog <- log(scale)
        for (i in 1:ns){
            haz[i, ] <- hlnorm(xx, meanlog = meanlog[i],
                               sdlog = sdlog[i]) * score
            sur[i, ] <- plnorm(xx, meanlog = meanlog[i],
                               sdlog = sdlog[i], lower.tail = FALSE)^score
            Haz[i, ] <- Hlnorm(xx, meanlog = meanlog[i],
                               sdlog = sdlog[i]) * score
        }
        dist = "Lognormal"
    }else if (x$dist == "ev"){
        for (i in 1:ns){
            haz[i, ] <- hEV(xx, shape = shape[i],
                                 scale = scale[i]) * score
            sur[i, ] <- pEV(xx, shape = shape[i],
                                 scale = scale[i], lower.tail = FALSE)^score
            Haz[i, ] <- HEV(xx, shape = shape[i],
                                 scale = scale[i]) * score
        }

        dist = "Extreme value"
    }else if (x$dist == "gompertz"){
        if (x$param == "canonical"){
            for (i in 1:ns){
                haz[i, ] <- hgompertz(xx, shape = score * shape[i],
                                      scale = scale[i],
                                      param = "canonical")## * score
                sur[i, ] <- pgompertz(xx, shape = score * shape[i],
                                      scale = scale[i],
                                      lower.tail = FALSE,
                                      param = "canonical") ##^score
                Haz[i, ] <- Hgompertz(xx, shape = score * shape[i],
                                      scale = scale[i],
                                      param = "canonical") ##* score
            }
        }else if (x$param == "rate"){
            for (i in 1:ns){
                haz[i, ] <- exp(shape + xx * scale) * score
                Haz[i, ] <- exp(shape) * score * expm1(xx * scale) / scale
                sur[i, ] <- exp(-Haz[, i])
            }
        }

        dist = "Gompertz"
    }

    if ("haz" %in% fn){

        if (is.null(ylim)) {
            ylim0 <- c(0, max(haz))
        }else{
            ylim0 <- ylim
        }
        ##if (min(p) < 1) ylim0[2] <- min(ylim0[2], max(haz[, -1]))

        if (is.null(xlab)) xlab <- "Duration"
        if (is.null(ylab)) ylab <- "Hazard"
        if (is.null(main)){
            hmain <- paste(dist, "hazard function")
        }else{
            hmain <- main
        }
        plot(xx, haz[1, ], type = "l", xlim = xlim, ylim = ylim0,
             col = col[1], lty = lty[1],
             xlab = xlab, ylab = ylab, main = hmain, ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, haz[i, ], type = "l", lty = lty[i], col = col[i])
            }
        }
        ##abline(h = 0)
        ##abline(v = 0)
        if ((ns > 1) && printLegend){
            legend(x = "bottomright",  legend = x$strata, lty = lty,
                   inset = 0.001,
                   col = col)
        }
    }
    ## Cumulative hazard
    if ("cum" %in% fn){

        if (is.null(ylim)){
            ylim0 <- c(0, max(Haz))
        }else{
            ylim0 <- ylim
            ylim0[2] <- max(ylim0[2], max(Haz))
        }
        ##if (is.null(xlab))
        xlab <- "Duration"
        ##if (is.null(ylab))
        ylab <- "Cumulative Hazard"
        if (is.null(main)){
            Hmain <- paste(dist, "cumulative hazard function")
        }else{
            Hmain <- main
        }
        plot(xx, Haz[1, ], type = "l", xlim = xlim, ylim = ylim0,
             xlab = xlab, ylab = ylab, main = Hmain, col = col[1],
             lty = lty[1], ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, Haz[i, ], type = "l", lty = lty[i], col = col[i])
            }
        }
        ##abline(h = 0)
        ##abline(v = 0)
        if ((ns > 1) && printLegend){
            legend(x = "topleft",  legend = x$strata, lty = lty,
                   col = col, inset = 0.001)
        }
    }
    ## density
    if ("den" %in% fn){

        den <- haz * sur
        ##if (is.null(ylim))
        ylim <- c(0, max(den))

        ##if (min(p) < 1) ylim[2] <- min(max(den[, -1]))

        ##if (is.null(xlab))
        xlab <- "Duration"
        ##if (is.null(ylab))
        ylab <- "Density"
        if (is.null(main)){
            dmain <- paste(dist, "density function")
        }else{
            dmain <- main
        }
        plot(xx, den[1, ], type = "l", xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, main = dmain, lty = lty[1],
             col = col[1], ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, den[i, ], type = "l", lty = lty[i], col = col[i])
            }
        }
        abline(h = 0)
        abline(v = 0)

        if ((ns > 1) && printLegend){
            legend(x = "topright",  legend = x$strata, lty = lty,
                   inset = 0.001,
                   col = col)
        }
        
    }
    ## Survivor function
    if ("sur" %in% fn){
        
        ##if (is.null(ylim))
        ylim <- c(0, 1)

        ##if (is.null(xlab))
        xlab <- "Duration"
        ##if (is.null(ylab))
        ylab <- "Survival"
        if (is.null(main)){
            smain <- paste(dist, "survivor function")
        }else{
            smain <- main
        }
        plot(xx, sur[1, ], type = "l", xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, main = smain, lty = lty[1],
             col = col[1], ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, sur[i, ], type = "l", lty = lty[i],
                      col = col[i])
            }
        }
        abline(h = 0)
        abline(v = 0)

        if ((ns > 1) && printLegend){
            legend(x = "bottomleft",  legend = x$strata, lty = lty,
                   inset = 0.001,
                   col = col)
        }
        

    }
    ##par(oldpar)
}
