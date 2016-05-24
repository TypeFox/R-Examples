"plot.bic.glm" <-
function (x, e = 1e-04, mfrow = NULL, include = 1:length(x$namesx), ...) 
{
    plotvar <- function(prob0, mixprobs, means, sds, Emean, Esd, 
        name, e = 1e-04, nsteps = 500, ...) {
        if (prob0 == 1) {
            xlower <- -0
            xupper <- 0
            xmax <- 1
        }
        else {
            qmin <- qnorm(e/2, Emean, Esd)
            qmax = qnorm(1 - e/2, Emean, Esd)
            xlower <- min(qmin, 0)
            xupper <- max(0, qmax)
        }
        xx <- seq(xlower, xupper, length.out = nsteps)
        yy <- rep(0, times = length(xx))
        maxyy <- 1
        if (prob0 < 1) {
            for (j in 1:length(means)) yy <- yy + mixprobs[j] * 
                dnorm(xx, means[j], sds[j])
            maxyy <- max(yy)
        }
        ymax <- max(prob0, 1 - prob0)
        plot(c(xlower, xupper), c(0, ymax), type = "n", xlab = "", 
            ylab = "", main = name, ...)
        lines(c(0, 0), c(0, prob0), lty = 1, lwd = 3, ...)
        lines(xx, (1 - prob0) * yy/maxyy, lty = 1, lwd = 1, ...)
    }
    vars <- unlist(x$assign[include + 1])
    nvar <- length(vars)
    probs <- NULL
    for (i in include) probs <- c(probs, rep(x$probne0[i], times = length(x$assign[[i + 
        1]])))
    nms <- NULL
    for (i in include) {
        if (is.na(x$output.names[i][1])) 
            nms <- c(nms, names(x$output.names[i]))
        else nms <- c(nms, paste(names(x$output.names[i]), unlist(x$output.names[i])[-1], 
            sep = "."))
    }
    wwhich <- NULL
    for (i in include) {
        wwhich <- cbind(wwhich, matrix(rep(x$which[, i], times = length(x$assign[[i + 
            1]])), ncol = length(x$assign[[i + 1]])))
    }
    if (!is.null(mfrow)) {
        lo <- mfrow
        losize <- lo[1] * lo[2]
    }
    else {
        layoutsizes <- c(1, 4, 9)
        layouts <- rbind(c(1, 1), c(2, 2), c(3, 3))
        layout <- max((1:length(layoutsizes))[layoutsizes <= 
            nvar])
        losize <- layoutsizes[layout]
        lo <- layouts[layout, ]
    }

    keep.mfrow = par()$mfrow
    par(mfrow = c(1, 1))
    par(ask = FALSE)
    par(mfrow = lo)
    ngroups <- ceiling(nvar/losize)
    for (k in 1:ngroups) {
        for (ii in ((k - 1) * losize + 1):min(k * losize, nvar)) {
            i <- vars[ii] - 1
            prob0 <- 1 - probs[i]/100
            sel <- wwhich[, i]
            mixprobs <- x$postprob[sel]/(1 - prob0)
            means <- x$mle[sel, i + 1]
            sds <- x$se[sel, i + 1]
            Emean <- x$condpostmean[i + 1]
            Esd <- x$condpostsd[i + 1]
            name <- nms[i]
            plotvar(prob0, mixprobs, means, sds, Emean, Esd, 
                name, e = e, ...)
        }
        par(ask = TRUE)
    }
   par(mfrow=keep.mfrow)
   par(ask = FALSE)
}
"plot.bicreg" <-
function (x, e = 1e-04, mfrow = NULL, include = 1:x$n.vars, include.intercept = TRUE, ...) 
{
    plotvar <- function(prob0, mixprobs, means, sds, Emean, Esd, 
        name, e = 1e-04, nsteps = 500, ...) {
        if (prob0 == 1) {
            xlower <- -0
            xupper <- 0
            xmax <- 1
        }
        else {
            qmin <- qnorm(e/2, Emean, Esd)
            qmax = qnorm(1 - e/2, Emean, Esd)
            xlower <- min(qmin, 0)
            xupper <- max(0, qmax)
        }
        xx <- seq(xlower, xupper, length.out = nsteps)
        yy <- rep(0, times = length(xx))
        maxyy <- 1
        if (prob0 < 1) {
            for (j in 1:length(means)) yy <- yy + mixprobs[j] * 
                dnorm(xx, means[j], sds[j])
            maxyy <- max(yy)
        }
        ymax <- max(prob0, 1 - prob0)
        plot(c(xlower, xupper), c(0, ymax), type = "n", xlab = "", 
            ylab = "", main = name, ...)
        lines(c(0, 0), c(0, prob0), lty = 1, lwd = 3, ...)
        lines(xx, (1 - prob0) * yy/maxyy, lty = 1, lwd = 1, ...)
    }
    nvar <- length(include) + include.intercept
    vars <- include + 1
    if (include.intercept) 
        vars <- c(1, vars)
    if (!is.null(mfrow)) {
        lo <- mfrow
        losize <- lo[1] * lo[2]
    }
    else {
        layoutsizes <- c(1, 4, 9)
        layouts <- rbind(c(1, 1), c(2, 2), c(3, 3))
        layout <- max((1:length(layoutsizes))[layoutsizes <= 
            nvar])
        losize <- layoutsizes[layout]
        lo <- layouts[layout, ]
    }

    keep.mfrow = par()$mfrow
    par(mfrow = c(1, 1))
    par(ask = FALSE)
    par(mfrow = lo)
    ngroups <- ceiling(nvar/losize)
    for (k in 1:ngroups) {
        for (ii in ((k - 1) * losize + 1):min(k * losize, nvar)) {
            i <- vars[ii]
            if (i != 1) {
                prob0 <- 1 - x$probne0[i - 1]/100
                sel <- x$which[, i - 1]
            }
            else {
                prob0 <- 0
                sel <- rep(TRUE, times = x$n.models)
            }
            mixprobs <- x$postprob[sel]/(1 - prob0)
            means <- x$ols[sel, i]
            sds <- x$se[sel, i]
            Emean <- x$condpostmean[i]
            Esd <- x$condpostsd[i]
            if (i == 1) 
                name <- "Intercept"
            else name <- x$namesx[i - 1]
            plotvar(prob0, mixprobs, means, sds, Emean, Esd, 
                name, e = e, ...)
        }
        par(ask = TRUE)
    }

   par(mfrow=keep.mfrow)
   par(ask = FALSE)
}
"plot.bic.surv" <-
function (x, e = 1e-04, mfrow = NULL, include = 1:length(x$namesx), ...) 
{
    plotvar <- function(prob0, mixprobs, means, sds, Emean, Esd, 
        name, e = 1e-04, nsteps = 500, ...) {
        if (prob0 == 1) {
            xlower <- -0
            xupper <- 0
            xmax <- 1
        }
        else {
            qmin <- qnorm(e/2, Emean, Esd)
            qmax = qnorm(1 - e/2, Emean, Esd)
            xlower <- min(qmin, 0)
            xupper <- max(0, qmax)
        }
        xx <- seq(xlower, xupper, length.out = nsteps)
        yy <- rep(0, times = length(xx))
        maxyy <- 1
        if (prob0 < 1) {
            for (j in 1:length(means)) yy <- yy + mixprobs[j] * 
                dnorm(xx, means[j], sds[j])
            maxyy <- max(yy)
        }
        ymax <- max(prob0, 1 - prob0)
        plot(c(xlower, xupper), c(0, ymax), type = "n", xlab = "", 
            ylab = "", main = name, ...)
        lines(c(0, 0), c(0, prob0), lty = 1, lwd = 3, ...)
        lines(xx, (1 - prob0) * yy/maxyy, lty = 1, lwd = 1, ...)
    }
    vars <- unlist(x$assign[include + 1])
    nvar <- length(vars)
    probs <- NULL
    for (i in include) probs <- c(probs, rep(x$probne0[i], times = length(x$assign[[i + 
        1]])))
    nms <- NULL
    for (i in include) {
        if (is.na(x$output.names[i][1])) 
            nms <- c(nms, names(x$output.names[i]))
        else nms <- c(nms, paste(names(x$output.names[i]), unlist(x$output.names[i])[-1], 
            sep = "."))
    }
    wwhich <- NULL
    for (i in include) {
        wwhich <- cbind(wwhich, matrix(rep(x$which[, i], times = length(x$assign[[i + 
            1]])), ncol = length(x$assign[[i + 1]])))
    }
    if (!is.null(mfrow)) {
        lo <- mfrow
        losize <- lo[1] * lo[2]
    }
    else {
        layoutsizes <- c(1, 4, 9)
        layouts <- rbind(c(1, 1), c(2, 2), c(3, 3))
        layout <- max((1:length(layoutsizes))[layoutsizes <= 
            nvar])
        losize <- layoutsizes[layout]
        lo <- layouts[layout, ]
    }

    keep.mfrow = par()$mfrow
    par(mfrow = c(1, 1))
    par(ask = FALSE)
    par(mfrow = lo)
    ngroups <- ceiling(nvar/losize)
    for (k in 1:ngroups) {
        for (ii in ((k - 1) * losize + 1):min(k * losize, nvar)) {
            i <- vars[ii]
            prob0 <- 1 - probs[i]/100
            sel <- wwhich[, i]
            mixprobs <- x$postprob[sel]/(1 - prob0)
            means <- x$mle[sel, i]
            sds <- x$se[sel, i]
            Emean <- x$condpostmean[i]
            Esd <- x$condpostsd[i]
            name <- nms[i]
            plotvar(prob0, mixprobs, means, sds, Emean, Esd, 
                name, e = e, ...)
        }
        par(ask = TRUE)
    }
   par(mfrow=keep.mfrow)
   par(ask = FALSE)
}


