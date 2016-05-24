## 2010-12-19 added argument noccasions in calls to pdot
## 2012-11-03 blocked gpclib
## 2012-12-18 warning for usage
## 2012-12-24 binomN argument
## 2013-04-19 hazard halfnormal, hazard exponential
## 2013-04-20 hazard annular normal, hazard cumulative gamma\
## 2013-04-23 hazard hazard rate
## 2014-02-13 removed gpclib
## 2014-03-12 bufferbiascheck() shifted from secr.fit
## 2014-08-25 comment out lth which was unused and called undefined fn get.pts

bias.D <- function (buffer, traps, detectfn, detectpar, noccasions, binomN = NULL, control = NULL) {
    gr <- function (r) {
        if (detectfn == 7) {
            CV2 <- (detectpar$z / detectpar$sigma)^2
            sdlog <- log(1 + CV2)^0.5
            meanlog <- log(detectpar$sigma) - sdlog^2/2
        }
        if (detectfn == 11) {
            mu <- detectpar$beta0 - 10 * log ( r^2 ) / 2.302585 + detectpar$beta1 * (r-1)
            mu[r<1] <- detectpar$beta0
        }

        switch (detectfn+1,
                detectpar$g0 * exp(-r^2/2/detectpar$sigma^2),
                detectpar$g0 * ( 1 - exp(-(r/detectpar$sigma)^-detectpar$z)),
                detectpar$g0 * exp(-r/detectpar$sigma),
                detectpar$g0 * (1 - (1 - exp(-r^2 / 2 / detectpar$sigma^2))^detectpar$z),
                detectpar$g0 * (r <= detectpar$sigma),
                detectpar$g0 * ( ifelse(r <= detectpar$w, rep(1,length(r)),
                                        exp(- (r - detectpar$w) / detectpar$sigma))),
                detectpar$g0 * exp(-(r-detectpar$w)^2/2/detectpar$sigma^2),
                detectpar$g0 * ( plnorm(r, meanlog, sdlog, lower.tail = FALSE)),
                detectpar$g0 * ( 1 - pgamma(r, shape = detectpar$z,
                    scale = detectpar$sigma/detectpar$z)),
                pnorm (-(detectpar$b0 + detectpar$b1 * r), mean=0, sd=1, lower.tail = FALSE),
                pnorm (q = detectpar$cutval, mean = detectpar$beta0 + detectpar$beta1 * r,
                   sd = detectpar$sdS, lower.tail = FALSE),
                pnorm (q = detectpar$cutval, mean = mu, sd = detectpar$sdS, lower.tail = FALSE),
                ,,
                1 - exp(-detectpar$lambda0 * exp(-r^2/2/detectpar$sigma^2)),
                1 - exp(-detectpar$lambda0 * (1 - exp(-(r/detectpar$sigma)^-detectpar$z))),
                1 - exp(-detectpar$lambda0 * exp(-r/detectpar$sigma)),
                1 - exp(-detectpar$lambda0 * exp(-(r-detectpar$w)^2/2/detectpar$sigma^2)),
                1 - exp(-detectpar$lambda0 * ( 1 - pgamma(r, shape = detectpar$z,
                    scale = detectpar$sigma/detectpar$z)))
          )
    }

    integrand1 <- function (r) {
        ## r < trapspacing/2
        if (control$method == 1)
            (1 - (1 - gr(r))^ (noccasions*k) )   * ntraps * 2 * pi *r
        else
            invlogit(predict(pdotr.spline,r)$y) * ntraps * 2 * pi *r
    }
    integrand2 <- function (r) {
        ## r >= trapspacing/2
        if (control$method == 1)
            (1 - (1 - gr(r))^ (noccasions*k) )   * l(r)
        else
            invlogit(predict(pdotr.spline,r)$y) * l(r)
    }

#    lth <- function (x) {
#        hull <- get.pts(x)
#        sum(sapply(hull, function(xy) {
#            xy$x <- c(xy$x, xy$x[1])
#            xy$y <- c(xy$y, xy$y[1])
#            sum(sqrt(diff(xy$x)^2 + diff(xy$y)^2))
#        }))
#    }
#    perimeterfn <- function (buffer, traps, ntheta = 60) {
#        theta <- (2 * pi) * (0:ntheta)/ntheta
#        pts <- data.frame(x = buffer * cos(theta), y = buffer * sin(theta))
#        centres <- split(traps,rownames(traps))
#        polys <- lapply(centres, function (x) {
#            temp <- sweep(pts, MARGIN = 2, FUN='+', STATS = unlist(x))
#        })
#        lth(punion2(polys, ntheta))
#    }

    detectfn <- valid.detectfn(detectfn)
    if (!(detector(traps) %in% .localstuff$pointdetectors))
        stop ("bias.D() requires passive point detectors (not polygon or transect)")
    if (!(detector(traps) %in% .localstuff$individualdetectors))
        stop ("bias.D() requires passive individual detectors (not unmarked or presence)")
    if (!is.null(usage(traps))) {
        if (any(usage(traps) != 1))
        warning ("bias.D() does not allow for variable effort (detector usage)")
    }

    ntraps <- nrow(traps)
    trapspacing <- spacing(traps)
    defaultcontrol <- list(bfactor = 20, masksample = 1000, spline.df = 10,
                           scale = 10000, ntheta = 60, method = 1)
    if (detectfn %in% c(1,2,7))
       defaultcontrol$bfactor <- 200
    if (is.null(control))
        control <- defaultcontrol
    else
        control <- replace (defaultcontrol, names(control), control)

    if (FALSE) {
        ## sidelined code
        buffs <- (round(bfactor):round(bfactor*2)) * trapspacing
        wayout <- sweep(cbind(buffs,rep(0,length(buffs))), STAT = unlist(apply(traps,2,max)),
                        MAR=2, FUN='+')
        pdotwo <- pdot(wayout, traps, detectfn, detectpar, noccasions, binomN)
        tempesa <- esa.plot(traps, max.buffer= trapspacing*control$bfactor,
                       spacing = trapspacing/2, detectfn = detectfn,
                       detectpar=detectpar, noccasions=noccasions,
                       plt = F, thin = 0.01)
        pdotr.spline <- smooth.spline(c(tempesa$buffer,buffs), c(tempesa$pdot, pdotwo),
                                      control$spline.df)
    }
    else {

## wasteful..
        temp <- make.mask(traps, buffer = trapspacing * control$bfactor,
                          spacing = trapspacing/2)
        OK <- sample(1:nrow(temp), size = min(nrow(temp), control$masksample))
        temp2 <- temp[OK,]

## new
#        temp2 <- matrix(runif (2 * control$masksample),nc=2) - 0.5
#        buff <- trapspacing * control$bfactor
#        dx <- diff(range(traps$x))
#        dy <- diff(range(traps$y))
#        temp2[,1] <- temp2[,1] * (dx + 2 * buff) + mean(traps$x)
#        temp2[,2] <- temp2[,2] * (dy + 2 * buff) + mean(traps$y)

        buff <- distancetotrap(temp2, traps)
        temp3 <- pdot(temp2, traps, detectfn, detectpar, noccasions, binomN)
        if (control$method == 1) {
            tempfit <- nls ( temp3 ~ (1 - (1 - gr(buff))^ (noccasions*k) ), start=list(k=2))
            if (tempfit$convInfo$isConv)
               k <- coef(tempfit)
            else
                stop ("failed to fit pdot vs buffer curve")
        }
        else if (control$method == 2) {
            temp3 <- logit(1 - (1 - temp3)^noccasions)
            OK <- is.finite(temp3)
            pdotr.spline <- smooth.spline(buff[OK],temp3[OK], df = control$spline.df)
        }
        else
            stop ("unrecognised 'control$method'")
    }

    ## make function to return linear approximation to contour length at radius r
    ## assuming spacing/2 <= r < (scale*spacing)

    hull <- buffer.contour(traps, buffer = 0, convex = TRUE, plt = FALSE, ntheta = 1)
    perimeter <- sum(sapply(hull, function (xy) sum(sqrt(diff(xy$x)^2 + diff(xy$y)^2))))
    perimeter <- perimeter + 2 * pi * trapspacing/2^0.5
    critx <- c(trapspacing/2,
               trapspacing/2^0.5,
               trapspacing * control$scale)  ## coarse
    crity <- c(ntraps * 2 * pi * trapspacing/2,
               perimeter,
               perimeter + 2 * pi * (trapspacing * control$scale - trapspacing/2^0.5))
    l <- approxfun(critx, crity, rule = 2)

    # scalar for 0-Inf
    I1 <- integrate (integrand1, lower = 0, upper = trapspacing/2)$value
    I2 <- integrate (integrand2, lower = trapspacing/2, upper = trapspacing/2^0.5)$value
    I3 <- integrate (integrand2, lower = trapspacing/2^0.5, upper = Inf)$value

    # in case buffer is vector
    nw <- length(buffer)
    I1w <- sapply(buffer, function(x) integrate (integrand1, lower = 0,
        upper = min(x,trapspacing/2))$value)
    I2w <- ifelse(buffer>trapspacing/2, sapply(buffer, function(x)
            integrate (integrand2, lower = trapspacing/2,
                upper =  min(x,trapspacing/2^0.5))$value), rep(0,nw))
    I3w <- ifelse(buffer>trapspacing/2^0.5, sapply(buffer, function(x)
            integrate (integrand2, lower = trapspacing/2^0.5,
                upper = x)$value), rep(0,nw))

    data.frame(buffer = buffer, RB.D = (I1+I2+I3)/ (I1w+I2w+I3w) - 1)
}

suggest.buffer <- function (object, detectfn = NULL, detectpar = NULL, noccasions = NULL,
    ignoreusage = FALSE,  RBtarget = 0.001, interval = NULL, binomN = NULL, ...) {
    if (ms(object)) {
        if (inherits(object,'secr')) {
            nsess <- length(object$capthist)
            traps <- traps(object$capthist)
            detectpar <- detectpar(object)
            detectfn <- object$detectfn
            noccasions <- sapply(object$capthist, ncol)
        }
        else {
            nsess <- length(object)
            if (inherits(object,'traps'))
                traps <- object
            if (is.null(detectpar)) {
                detectpar <- vector(nsess, mode='list')
            }
            else {
                if (length(detectpar[[1]]) == 1)
                    detectpar <- rep(detectpar, nsess)
            }
            if (length(noccasions) == 1)
                noccasions <- rep(noccasions, nsess)
        }

        buffer <- numeric(nsess)
        for (i in 1:nsess) {
            if (inherits(object,'capthist'))
                buffer[i] <- suggest.buffer(object[[i]], detectfn,
                    detectpar[[i]], noccasions[i], ignoreusage, RBtarget, interval, binomN, ...)
            else
                buffer[i] <- suggest.buffer(traps[[i]], detectfn,
                    detectpar[[i]], noccasions[i], ignoreusage, RBtarget, interval, binomN, ...)
         }
        buffer
    }
    else {
        if (inherits(object,'secr')) {
            traps <- traps(object$capthist)
            if (!is.null(detectfn))
                if (detectfn != object$detectfn)
                    warning ("changing 'detectfn' is risky", call.=FALSE)
            if (is.null(detectpar))
                detectpar <- detectpar(object)
            if (is.null(detectfn))
                detectfn <- object$detectfn
            if (is.null(noccasions))
                noccasions <- ncol(object$capthist)
            if (is.null(ignoreusage)) {
                if (!is.null(object$details$ignoreusage))
                    ignoreusage <- object$details$ignoreusage   ## assume not an old model
                else
                    ignoreusage <- FALSE
            }
        }
        else {
            if (inherits(object, 'capthist')) {
                traps <- traps(object)
                noccasions <- ncol(object)
                if (is.null(detectpar)) {
                    tempmask <- make.mask (traps, 6*RPSV(object))
                    detectpar <- autoini (object, tempmask,
                        ignoreusage = ignoreusage)[parnames(0)]
                    tempdp <- lapply(detectpar,formatC,4)
                    warning ("using automatic 'detectpar' ",
                        paste(names(detectpar), "=", tempdp, collapse=", "),
                        call. = FALSE)
                    if (!is.null(detectfn))
                        if (detectfn != 0)
                            warning ("forcing 'detectfn' to halfnormal", call. = FALSE)
                    detectfn <- 0
                }
            }
            else {
                traps <- object
                ## could retrieve noccasions from usage here...
            }
        }
        if (is.null(interval)) {
            interval <- c(1, 100 * spatialscale(detectpar, detectfn))
        }
        detectfn <- valid.detectfn(detectfn)
        detectpar <- valid.detectpar(detectpar, detectfn)
        if (!(detector(traps) %in% .localstuff$pointdetectors))
            stop ("require passive point detectors (not polygon or transect)")
        if (!(detector(traps) %in% .localstuff$individualdetectors))
            stop ("require passive individual detectors (not unmarked or presence)")

        if (ignoreusage)
            usage(traps) <- NULL
        fn <- function (w) {
            bias.D(w, traps, detectfn, detectpar, noccasions, binomN, ...)$RB.D - RBtarget
        }
        temp <- try(round(uniroot (fn, interval)$root), silent = TRUE)
        if (inherits(temp, 'try-error')) {
            stop("buffer outside interval ",
                 formatC(interval[1],5),' to ', formatC(interval[2],5), " m")
        }
        else
            temp
    }
}

## shifted from secr.fit 2014-03-12

bufferbiascheck <- function (output, buffer, biasLimit) {

    ############################################
    ## buffer bias check
    ## not for polygon & transect detectors
    ############################################

    capthist <- output$capthist
    validbiasLimit <- !is.null(biasLimit)
    validbiasLimit <- validbiasLimit & is.finite(biasLimit)
    validbiasLimit <- validbiasLimit & (biasLimit>0)
    if ((output$fit$value < 1e9) &
        (detector(traps(capthist)) %in% .localstuff$pointdetectors) &
        !(detector(traps(capthist)) %in% c('unmarked','presence')) &
        is.null(telemetryxy(capthist)) &
        validbiasLimit) {
        if (ms(capthist)) {
            nsess <- length(capthist)
            bias <- numeric(nsess)
            for (i in 1:nsess) {
                temptrps <- traps(capthist)[[i]]
                if (output$details$ignoreusage)
                    usage(temptrps) <- NULL
                dpar <-  detectpar(output)[[i]]
                biastemp <- try( bias.D(buffer, temptrps,
                                        detectfn = output$detectfn,
                                        detectpar = dpar,
                                        noccasions = ncol(capthist[[i]]),
                                        binomN = output$details$binomN) )
                if (inherits(biastemp, 'try-error'))
                    warning('could not perform bias check')
                else
                    bias[i] <- biastemp$RB.D
            }
        }
        else {
            temptrps <- traps(capthist)
            if (output$details$ignoreusage)
                usage(temptrps) <- NULL
            dpar <-  detectpar(output)
            bias <- try( bias.D(buffer, temptrps,
                                detectfn = output$detectfn,
                                detectpar = dpar,
                                noccasions = ncol(capthist),
                                binomN = output$details$binomN) )
            if (inherits(bias, 'try-error')) {
                warning('could not perform bias check')
                bias <- 0  ## 2012-12-18 suppresses second message
            }
            else
                bias <- bias$RB.D
        }
        if (any(bias > biasLimit))
            warning ("predicted relative bias exceeds ", biasLimit, " with ",
                     "buffer = ", buffer)
    }
}
