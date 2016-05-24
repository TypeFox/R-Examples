################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Plot estimated interaction kernel (siaf/tiaf) as a function of distance
###
### Copyright (C) 2012-2015 Sebastian Meyer
### $Revision: 1325 $
### $Date: 2015-04-30 13:56:23 +0200 (Don, 30. Apr 2015) $
################################################################################


iafplot <- function (object, which = c("siaf", "tiaf"), types = NULL,
    scaled = c("intercept", "standardized", "no"), truncated = FALSE, log = "",
    conf.type = if (length(pars) > 1) "MC" else "parbounds",
    conf.level = 0.95, conf.B = 999,
    xgrid = 101, col.estimate = rainbow(length(types)), col.conf = col.estimate,
    alpha.B = 0.15, lwd = c(3,1), lty = c(1,2),
    verticals = FALSE, do.points = FALSE,
    add = FALSE, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL,
    legend = !add && (length(types) > 1), ...)
{
    if (isTRUE(verticals)) verticals <- list()
    if (isTRUE(do.points)) do.points <- list()
    if (add) log <- paste0("", if (par("xlog")) "x", if (par("ylog")) "y")
    scaled <- if (is.logical(scaled)) { # surveillance < 1.9-0
        if (scaled) "intercept" else "no"
    } else {
        match.arg(scaled)
    }
    coefs <- coef(object)
    epiloglink <- .epilink(object) == "log"
    typeNames <- rownames(object$qmatrix)
    nTypes <- length(typeNames)
    
    ## interaction function
    which <- match.arg(which)
    IAFobj <- object$formula[[which]]
    if (is.null(IAFobj))
        stop("the model has no epidemic component")
    IAF <- IAFobj[[if (which=="siaf") "f" else "g"]]
    if (which == "siaf") { # needs to be a function of distance
        IAF <- as.function(
            c(alist(x=, ...=), quote(f(cbind(x, 0), ...))),
            envir = list2env(list(f = IAF), parent = environment(IAF))
            )
    }
    isStepFun <- !is.null(knots <- attr(IAFobj, "knots")) &&
        !is.null(maxRange <- attr(IAFobj, "maxRange"))

    ## interaction range
    if (isScalar(truncated)) {
        eps <- truncated
        truncated <- TRUE
    } else {
        eps <- attr(IAFobj, "eps")
    }
    if (is.null(eps)) { # cannot take eps into account (pre 1.8-0 behaviour)
        eps <- NA_real_
    } else if (length(eps) > 1L && truncated) {
        message("no truncation due to heterogeneous interaction ranges, see \"rug\"")
    }
    epsIsFixed <- length(eps) == 1L && is.finite(eps)

    ## scaled interaction function
    if (scaled == "intercept") {
        idxgamma0 <- match("e.(Intercept)", names(coefs), nomatch = 0L)
        if (idxgamma0 == 0L) {
            message("no scaling due to missing epidemic intercept")
            scaled <- "no"
        }
    } else { # we do not use gamma0 -> 0-length selection
        idxgamma0 <- 0L
    }
    SCALE <- switch(scaled,
        "intercept" = if (epiloglink) quote(exp(gamma0)) else quote(gamma0),
        "standardized" = quote(1/IAF(0, iafpars, types)),
        "no" = 1
    )
    FUN <- function (x, iafpars, types, gamma0) {
        scale <- eval(SCALE)
        vals <- scale * IAF(x, iafpars, types)
    }
    
    ## truncate at eps
    if (truncated && epsIsFixed) {
        body(FUN) <- as.call(c(as.list(body(FUN)), expression(
            vals[x > eps] <- 0,
            vals
            )))
    }
    
    ## if (loglog) {
    ##     body(FUN)[[length(body(FUN))]] <-
    ##         call("log", body(FUN)[[length(body(FUN))]])
    ## }

    ## extract parameters
    gamma0 <- coefs[idxgamma0]
    idxiafpars <- grep(paste0("^e\\.",which), names(coefs))
    iafpars <- coefs[idxiafpars]
    ## concatenate parameters
    idxpars <- c(idxgamma0, idxiafpars)
    pars <- c(gamma0, iafpars)

    ## type of confidence band
    force(conf.type)                    # per default depends on 'pars'
    if (length(pars) == 0 || is.null(conf.type) || is.na(conf.type)) {
        conf.type <- "none"
    }
    conf.type <- match.arg(conf.type,
                           choices = c("parbounds", "bootstrap", "MC", "none"))
    if (conf.type == "bootstrap") conf.type <- "MC"  # "bootstrap" was used <1.8
    if (conf.type == "parbounds" && length(pars) > 1) {
        warning("'conf.type=\"parbounds\"' is only valid for a single parameter")
    }

    ## grid of x-values (t or ||s||) on which FUN will be evaluated
    if (is.null(xlim)) {
        xmax <- if (add) {
            xmax <- par("usr")[2] / (if (par("xaxs")=="r") 1.04 else 1)
            if (par("xlog")) 10^xmax else xmax
        } else {
            if (epsIsFixed) {
                eps
            } else if (isStepFun && maxRange < Inf) {
                maxRange
            } else if (which == "siaf") {
                sqrt(sum((object$bbox[,"max"] - object$bbox[,"min"])^2))
            } else {
                diff(object$timeRange)
            }
        }
        xlim <- c(0.5*grepl("x", log), xmax)
    }
    xgrid <- if (isStepFun) {
        c(if (grepl("x", log)) {
            if (xlim[1L] < knots[1L]) xlim[1L] else NULL
        } else 0, knots[knots<xlim[2L]])
    } else if (isScalar(xgrid)) {
        if (grepl("x", log)) {
            exp(seq(log(xlim[1L]), log(xlim[2]), length.out=xgrid))
        } else seq(xlim[1L], xlim[2L], length.out=xgrid)
    } else {
        stopifnot(!is.na(xgrid), is.vector(xgrid, mode="numeric"))
        ## xgrid-specification overrides default xlim
        if (is.null(match.call()$xlim)) xlim <- range(xgrid)
        sort(xgrid)
    }

    ## type selection
    if (is.null(types)) { # check if the interaction function is type-specific
        if (nTypes == 1L || IAFobj$npars < 2) {
            types <- 1L
        } else { # we compare the values for different types
            .fByType <- sapply(seq_len(nTypes), function (type)
                               FUN(xgrid, iafpars, type, gamma0))
            types <- if (all(apply(.fByType[,-1L,drop=FALSE], 2L,
                                   function (fvals)
                                   identical(.fByType[,1L], fvals))))
                1L else seq_len(nTypes)
        }
    }

    ## initialize plotting frame
    if (!add) {
        if (is.null(ylim))
            ylim <- if (grepl("y", log)) {
                sort(FUN(xlim, iafpars, 1L, gamma0))
            } else c(0, FUN(0, iafpars, 1L, gamma0))
        if (is.null(xlab)) xlab <- if (which == "siaf") {
            expression("Distance " * x * " from host")
        } else {
            expression("Time " * t * " since infectious")
        }
        if (is.null(ylab)) {
            ylab <- if (which == "siaf") {
                expression(f(x))        # f(group("||",bold(s)-bold(s)[j],"||"))
            } else {
                expression(g(t))
            }
            ## if (loglog) ylab[[1]] <- substitute(log(ylab), list(ylab=ylab[[1]]))
            if (scaled == "intercept") {
                ylab <- ##if (loglog) as.expression(call("+", quote(gamma[0]), ylab[[1]])) else
                    as.expression(call("paste",
                        if (epiloglink) quote(e^{gamma[0]}) else quote(gamma[0]),
                        quote(phantom() %.% phantom()), ylab[[1]]))
            }
            if (scaled == "standardized") {
                ylab <- as.expression(call("/", ylab[[1]],
                    if (which == "siaf") quote(f(0)) else quote(g(0))))
            }
        }
        plot(xlim, ylim, type="n", xlab = xlab, ylab = ylab, log = log, ...)
        if (length(eps) > 1L && truncated) rug(eps)
    }

    ## store evaluated interaction function in a matrix (will be returned)
    typeNamesSel <- typeNames[types]
    res <- matrix(NA_real_, length(xgrid), 1L+length(types),
                  dimnames = list(NULL, c("x", typeNamesSel)))
    res[,1L] <- xgrid
    for (i in seq_along(types)) {
        ## select parameters on which to evaluate iaf
        parSample <- switch(conf.type, parbounds = {
            cis <- confint(object, idxpars, level=conf.level)
            ## all combinations of parameter bounds
            do.call("expand.grid", as.data.frame(t(cis)))
        }, MC = { # Monte-Carlo confidence interval
            ## sample parameters from their asymptotic multivariate normal dist.
            rbind(pars,
                  mvrnorm(conf.B, mu=pars,
                          Sigma=vcov(object)[idxpars,idxpars,drop=FALSE]),
                  deparse.level=0)
        })
        
        ## add confidence limits
        if (!is.null(parSample)) {
            fvalsSample <- apply(parSample, 1, if (scaled == "intercept") {
                function (pars) FUN(xgrid, pars[-1L], types[i], pars[1L])
            } else {
                function (pars) FUN(xgrid, pars, types[i])
            })
            if (length(xgrid) == 1L)  # e.g., single-step function
                fvalsSample <- t(fvalsSample)  # convert to matrix form
            lowerupper <- if (conf.type == "parbounds") {
                t(apply(fvalsSample, 1, range))
            } else { # Monte-Carlo sample of parameter values
                if (is.na(conf.level)) {
                    stopifnot(alpha.B >= 0, alpha.B <= 1)
                    .col <- col2rgb(col.conf[i], alpha=TRUE)[,1]
                    .col["alpha"] <- round(alpha.B*.col["alpha"])
                    .col <- do.call("rgb", args=c(as.list(.col),
                                           maxColorValue = 255))
                    matlines(x=xgrid, y=fvalsSample, type="l", lty=lty[2],
                             col=.col, lwd=lwd[2]) # returns NULL
                } else {
                    t(apply(fvalsSample, 1, quantile,
                            probs=c(0,conf.level) + (1-conf.level)/2))
                }
            }
            if (!is.null(lowerupper)) {
                attr(res, if(length(types)==1) "CI" else
                     paste0("CI.",typeNamesSel[i])) <- lowerupper
                if (isStepFun) {
                    segments(rep.int(xgrid,2L), lowerupper,
                             rep.int(c(xgrid[-1L], min(maxRange, xlim[2L])), 2L), lowerupper,
                             lty=lty[2], col=col.conf[i], lwd=lwd[2])
                    ##points(rep.int(xgrid,2L), lowerupper, pch=16, col=col.conf[i])
                } else {
                    matlines(x=xgrid, y=lowerupper,
                             type="l", lty=lty[2], col=col.conf[i], lwd=lwd[2])
                }
            }
        }
        
        ## add point estimate
        res[,1L+i] <- FUN(xgrid, iafpars, types[i], gamma0)
        if (isStepFun) {
            segments(xgrid, res[,1L+i],
                     c(xgrid[-1L], min(maxRange, xlim[2L])), res[,1L+i],
                     lty = lty[1], col = col.estimate[i], lwd = lwd[1])
            ## add points
            if (is.list(do.points)) {
                pointStyle <- modifyList(list(pch=16, col=col.estimate[i]),
                                         do.points)
                do.call("points", c(list(xgrid, res[,1L+i]), pointStyle))
            }
            ## add vertical connections:
            if (is.list(verticals)) {
                verticalStyle <- modifyList(
                    list(lty = 3, col = col.estimate[i], lwd = lwd[1L]),
                    verticals)
                do.call("segments", c(
                    list(xgrid[-1L], res[-length(xgrid),1L+i],
                         xgrid[-1L], res[-1L,1L+i]), verticalStyle))
            }
            if (maxRange <= xlim[2L]) { ## add horizontal=0 afterwards
                segments(maxRange, 0, xlim[2L], 0,
                         lty = lty[1], col = col.estimate[i], lwd = lwd[1])
                if (is.list(verticals))
                    do.call("segments", c(
                        list(maxRange, res[length(xgrid),1L+i],
                             maxRange, 0), verticalStyle))
                if (is.list(do.points))
                    do.call("points", c(list(maxRange, 0), pointStyle))
            }
        } else {
            lines(x = xgrid, y = res[,1L+i],
                  lty = lty[1], col = col.estimate[i], lwd = lwd[1])
        }
    }
    
    ## add legend
    if (isTRUE(legend) || is.list(legend)) {
        default.legend <- list(x = "topright", legend = typeNamesSel,
                               col = col.estimate, lty = lty[1], lwd = lwd[1],
                               bty = "n", cex = 0.9, title="type")
        legend.args <- if (is.list(legend)) {
            modifyList(default.legend, legend)
        } else default.legend
        do.call("legend", legend.args)
    }

    ## Invisibly return interaction function evaluated on xgrid (by type)
    invisible(res)
}
