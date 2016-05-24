############################################################################################
## plot.secr.R
## Method for plotting detection function from fitted secr object
## 2009 02 12 2009 08 09 2009 09 13 (limits)
## 2009 11 04 signal strength
## 2010 07 01 alpha detectfn
## 2010 09 15 amend message for missing cutval
## 2010 11 01 enabled plot detectfn=9
## 2011 03 26 reversed order of first two arguments of dfn's
## default of limits changed to FALSE in plot.secr()
## 2012-10-21 moved HN etc to utility.r
## 2013-04-20 new detection functions HHN, HHR, HEX, HAN, HCG
## 2013-04-24 uses getdfn from utility.r
## 2014-04-25 ylim bug when g0 fixed
## 2014-09-18 limits = TRUE now works with acoustic dfn
############################################################################################

plot.secrlist <- function (x, newdata=NULL, add = FALSE,
    sigmatick = FALSE, rgr = FALSE, limits = FALSE, alpha = 0.05, xval = 0:200,
    ylim = NULL, xlab = NULL, ylab = NULL, ..., overlay = TRUE)
{
    if (overlay) {
        plot(x[[1]], newdata, add, sigmatick, rgr, limits, alpha, xval, ylim, xlab,
               ylab, ...)
        lapply(x[-1], plot, newdata, TRUE, sigmatick, rgr, limits, alpha, xval, ylim, xlab,
               ylab, ...)
    }
    else {
        lapply(x, plot, newdata, TRUE, sigmatick, rgr, limits, alpha, xval, ylim, xlab,
               ylab, ...)
    }
    invisible()
}

plot.secr <- function (x, newdata=NULL, add = FALSE,
    sigmatick = FALSE, rgr = FALSE, limits = FALSE, alpha = 0.05, xval = 0:200,
    ylim = NULL, xlab = NULL, ylab = NULL, ...)
{
    gline <- function (predicted, rowi = 1, eps = 1e-10) {
        ## eps is used to limit y to range where gradient() works
        ## may need later adjustment
        if (!is.data.frame(predicted)) {
            out <- list()
            for (i in 1:length(predicted))
                out[[i]] <- gline(predicted[[i]], i)
            names(out) <- names(predicted)
            out
        }
        else {
            pars <- predicted[parnames(x$detectfn),'estimate']
            pars[is.na(pars)] <- unlist(x$fixed)
            dfn <- getdfn(x$detectfn)
            ## superceded 2013-04-24
            ## switch (x$detectfn+1, HN, HR, EX, CHN, UN, WEX, ANN, CLN, CG, BSS, SS,
            ##               SSS, SN, SNS, HHN, HHR, HEX, HAN, HCG)
            if (sigmatick) {
              sigma <- pars[2]
              y <- dfn(sigma, pars, x$details$cutval)
              dy <- par()$usr[4]/20
              segments (sigma, y-dy, sigma, y+dy)
            }

            y <- dfn(xval, pars, x$details$cutval)

            if (rgr) {
              y <- xval * y
              ymax <- par()$usr[4]
              lines (xval, y * 0.8 * ymax / max(y), lty = 2, ...)
            }
            else lines (xval, y, ...)

            if (limits & !rgr) {
                ## delta method variance of g()

                grad <- matrix(nrow = length(xval), ncol = length(x$fit$par))  ## beta pars
                if (is.null(newdata)) newdata <- secr.make.newdata (x)

                lkdfn <- function (beta, r) {
                    ## real g() from all beta pars and model.matrix
                    parnamvec <- parnames(x$detectfn)
                    real <- numeric(length(parnamvec))
                    names(real) <- parnamvec

                    for (rn in parnamvec) {
                         par.rn <- x$parindx[[rn]]
                         ## 2014-08-19
                         ## mat <- model.matrix(x$model[[rn]], data=newdata[rowi,,drop=F])
                         mat <- general.model.matrix(x$model[[rn]], data=newdata[rowi,,drop=F])
                         lp <- mat %*% matrix(beta[par.rn], ncol = 1)
                         real[rn] <- untransform (lp, x$link[[rn]])
                    }
                    ## bug fix 2014-09-18 : requires x$details$cutval
                    logit(dfn(r, real, x$details$cutval))
                }

                for (i in 1:length(xval))

            ## Fast special cases: checking only
            ##    if ((x$detectfn==0) & (all(sapply(x$model, function(m) m == ~1)) ))
            ##    {
            ##        ## ASSUME DEFAULT LINK
            ##        g0 <- logit(pars[1])
            ##        sigma <- log(pars[2])
            ##        r <- xval[i]
            ##        ## D(expression(1/(1+exp(-g0)) * exp(-r^2/2/exp(sigma)^2)), 'g0')
            ##        ## D(expression(1/(1+exp(-g0)) * exp(-r^2/2/exp(sigma)^2)), 'sigma')
            ##        tempgrad <- c(exp(-g0)/(1 + exp(-g0))^2 * exp(-r^2/2/exp(sigma)^2),
            ##            -(1/(1 + exp(-g0)) * (exp(-r^2/2/exp(sigma)^2) * (-r^2/2 * (2 * (exp(sigma)
            ##            * exp(sigma)))/(exp(sigma)^2)^2))))
            ##        if (!x$CL) tempgrad <- c(0,tempgrad)  ## for density
            ##        grad[i,] <- tempgrad
            ##    }
            ##    else
            ##       if ((x$detectfn==1) & (all(sapply(x$model, function(m) m == ~1)) ))
            ##    {
            ##        ## ASSUME DEFAULT LINK
            ##        g0 <- logit(pars[1])
            ##        sigma <- log(pars[2])
            ##        z <- log(pars[3])
            ##           r <- xval[i]
            ##        tempgrad <- c(
            ##            exp(-g0)/(1 + exp(-g0))^2 * (1 - exp(-(r/exp(sigma))^(-exp(z)))),
            ##            -(1/(1 + exp(-g0)) * (exp(-(r/exp(sigma))^(-exp(z))) * ((r/exp(sigma))^((-exp(z)) -
            ##            1) * ((-exp(z)) * (r * exp(sigma)/exp(sigma)^2))))),
            ##            -(1/(1 + exp(-g0)) * (exp(-(r/exp(sigma))^(-exp(z))) * ((r/exp(sigma))^(-exp(z)) *
            ##            (log((r/exp(sigma))) * exp(z)))))
            ##        )
            ##        if (!x$CL) tempgrad <- c(0,tempgrad)  ## for density
            ##        grad[i,] <- tempgrad
            ##    }
            ##    else

                # grad[i,] <- fdHess (pars = x$fit$par, fun = lkdfn, r = xval[i])$gradient
                # grad[i,] <- grad (func = lkdfn, x = x$fit$par, r = xval[i])  ## needs numDeriv
                grad[i,] <- gradient (pars = x$fit$par, fun = lkdfn, r = xval[i])  ## see 'utility.R'

                vc <- vcov (x)
                gfn <- function(gg) {
                    gg <- matrix(gg, nrow = 1)
                    gg %*% vc %*% t(gg)
                    }
                se <- apply(grad, 1, gfn)^0.5
                ## lcl <- pmax(y - z*se,0)  # on natural scale
                ## ucl <- pmin(y + z*se,1)

                ## limits on link scale
                lcl <- ifelse ((y>eps) & (y<(1-eps)), invlogit (logit(y) - z*se), NA)
                ucl <- ifelse ((y>eps) & (y<(1-eps)), invlogit (logit(y) + z*se), NA)

                lines (xval, lcl, lty=2, ...)
                lines (xval, ucl, lty=2, ...)
            }

            if (limits & !rgr)
                data.frame(x=xval, y=y, lcl = lcl, ucl = ucl)
            else
                data.frame(x=xval, y=y)
        }
    }
    z <- abs(qnorm(1-alpha/2))   ## beware confusion with hazard-rate z!
    temp <- predict (x, newdata)
    if (is.null(ylim)) {
        if (x$detectfn %in% c(9,10,11,12,13)) {      ## included 9 2010-11-01
            ylim <- c(0, 1)
        }
        else {
            if (x$detectfn %in% 14:18)
                yname <- 'lambda0'
            else
                yname <- 'g0'
            getmax <- function(x) {
                g0 <- x[yname,'estimate']
                se.g0 <- x[yname,'SE.estimate']
                if (limits & is.finite(se.g0))  ## is.finite 2010-10-10
                    min(1, g0 + z * se.g0)
                else
                    g0
            }
            if (is.data.frame(temp)) maxg0 <- getmax(temp)
            else maxg0 <- max(sapply (temp, getmax))

            if (is.na(maxg0)) maxg0 <- x$fixed[[yname]]
            if (maxg0 > 0.75) maxg0 <- max(maxg0,1)
            ylim <- c(0, maxg0)
        }
    }
    if (!add) {
        if (is.null(xlab))
            xlab <- 'Distance  (m)'
        if (is.null(ylab)) {
           binomN <- ifelse(is.null(x$details$binomN),0,x$details$binomN)
## revisit this 2011-02-06
           dlambda <- (detector(traps(x$capthist)) %in% c('polygon','polygonX')) |
               ((detector(traps(x$capthist)) %in% c('count')) & (binomN==0))
           if (dlambda)
               ylab <- 'Detection lambda'
           else
               ylab <- 'Detection probability'
        }
        plot (type ='n', 0,0, xlim=range(xval), ylim=ylim,
            xlab=xlab, ylab=ylab,...)
    }
    invisible(gline(temp))
}
############################################################################################

detectfnplot <- function (detectfn, pars, details = NULL,
    add = FALSE, sigmatick = FALSE, rgr = FALSE,
    xval = 0:200, ylim = NULL, xlab = NULL, ylab = NULL, ...)
{

    gline <- function (pars) {
        ## here pars is a vector of parameter values
        dfn <- getdfn(detectfn)
        if (sigmatick) {
            sigma <- pars[2]
            y <- dfn(sigma, pars,details$cutval)
            dy <- par()$usr[4]/20
            segments (sigma, y-dy, sigma, y+dy)
        }
        y <- dfn(xval, pars, details$cutval)
        if (rgr) {
            y <- xval * y
            ymax <- par()$usr[4]
            lines (xval, y * 0.8 * ymax / max(y), lty = 2, ...)
        }
        else lines (xval, y, ...)

        data.frame(x=xval, y=y)

    }

    ### mainline

    if (is.list(pars)) {   ## 2010-10-26
        if (is.list(pars[[1]]))
            pars <- matrix(unlist(pars), nrow = length(pars), byrow = T)
        else
            pars <- unlist(pars)
    }

    if (!is.matrix(pars)) pars <- matrix(pars, nrow = 1)

    ## added 2010-07-01
    if (is.character(detectfn))
        detectfn <- detectionfunctionnumber(detectfn)

    needp <- c(2,3,2,3,2,3,3,3,3,2,3,3,5,5,2,3,2,3,3)[detectfn+1]

    if (ncol(pars) != needp)
        stop ("require ", needp, " parameters for ",
             detectionfunctionname(detectfn), " detection function")

    if (is.null(ylim)) {
        if (detectfn %in% c(10,11,12,13)) {
            ylim <- c(0, 1)
        }
        else {
            ylim <- c(0, max(pars[,1]))   ## g0 or lambda0
        }
    }

    if (!add) {
        if (is.null(xlab))
            xlab <- 'Distance  (m)'
        if (is.null(ylab))
            ylab <- 'Detection'
        plot (type ='n', 0,0, xlim=range(xval), ylim=ylim,
            xlab=xlab, ylab=ylab,...)
    }

    invisible( apply (pars, 1, gline) )

}
############################################################################################

attenuationplot <- function (pars, add = FALSE, spherical = TRUE,
    xval = 0:200, ylim = NULL, xlab = NULL, ylab = NULL, ...) {

    mufn <- function (pars, r) {
        beta0 <- pars[1]
        beta1 <- pars[2]
        ## if spherical, assume distance r measured from 1 m
        if (spherical) {
            mu <- beta0 - 10 * log ( r^2 ) / 2.302585 + beta1 * (r-1)
            mu[r<1] <- beta0
        }
        else
            mu <- beta0 + beta1 * r
        mu
    }

    aline <- function (pars, r) {
        y <- mufn (pars, xval)
        lines (xval, y, ...)
        data.frame(x=xval, y=y)
    }

    if (is.list(pars)) {   ## 2010-10-26
        if (is.list(pars[[1]]))
            pars <- matrix(unlist(pars), nrow = length(pars), byrow = T)
        else
            pars <- unlist(pars)
    }
    if (!is.matrix(pars)) pars <- matrix(pars, nrow = 1)

    if (is.null(ylim)) {
        lower <- min(apply(pars, 1, mufn, r = max(xval)))
        ylim <- c(lower, max(pars[,1]))
    }

    if (!add) {
        if (is.null(xlab))
            xlab <- 'Distance  (m)'
        if (is.null(ylab))
            ylab <- 'Acoustic power (dB)'
        plot (type ='n', 0,0, xlim=range(xval), ylim=ylim,
            xlab=xlab, ylab=ylab,...)
    }
    invisible( apply (pars, 1, aline) )

}

