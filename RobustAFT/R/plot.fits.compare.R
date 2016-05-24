plot.fits.compare <- function (x, xplots = FALSE, ask = TRUE, which = 1:4, leg.position = c("topleft", "topleft", "topleft"), ...)
{
    if(!inherits(x, "fits.compare"))
        stop("Use only with 'TML' objects")
    show <- rep(FALSE, 4)
    if(!is.numeric(which) || any(which < 1) || any(which > 4))
        stop("which must be in 1:4")
    show[which] <- TRUE
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    oldmfcol <- par()$mfcol
    par(mfcol = c(1, 1))
    if(length(which)==1 & xplots==FALSE) ask <- FALSE
    op <- par(ask = ask)
    names <- rep("", length(x))
    form <- vector("list", length(x))
    f <- vector("list", length(x))
    r <- vector("list", length(x))
    rs <- vector("list",length(x))
    den <- vector("list", length(x))
    y <- vector("list", length(x))
    q <- vector("list", length(x))
    yname <- vector("list", length(x))
    fname <- vector("list", length(x))
    sigma <- vector("list", length(x))
    denmax <- 0
    denmin <- 1e+20
    if (xplots) {
        indep <- vector("list", length(x))
    }
    for (i in 1:length(x)) {
        xx <- class(x[[i]])[1]
        if (xx == "TML") {
            x[[i]]$coefficients <- x[[i]]$th1
        }
    }
    if (xplots) {
        aname <- function(x) {
            names(x$coefficients)
        }
        unames <- unique(as.vector(sapply(x, aname)))
        k <- match("(Intercept)", unames, nomatch = 0)
        if (k)
            unames <- unames[-k]
        lnames <- length(unames)
        xmin <- rep(Inf, lnames)
        xmax <- rep(-Inf, lnames)
    }
    for (i in 1:length(x)) {
        xx <- class(x[[i]])[1]
        if (xx == "glm" || xx == "cubinf")
            sigma[[i]] <- x[[i]]$dispersion
        else if (xx == "TML") {
            sigma[[i]] <- x[[i]]$v1
        }
        else {
            zs <- x[[i]]$scale
            if (is.null(zs))
                zs <- sqrt(sum(residuals(x[[i]])^2)/x[[i]]$df.residual)
            sigma[[i]] <- zs
        }
        names[i] <- x[[i]]$name
        n <- length(residuals(x[[i]])) + 0.5
        form[[i]] <- formula(x[[i]])
        f[[i]] <- fitted.values(x[[i]])
        r[[i]] <- residuals(x[[i]])
        rs[[i]] <- r[[i]] / sigma[[i]]
        den[[i]] <- density(residuals(x[[i]]))
        denmax <- max(denmax, den[[i]]$y)
        denmin <- min(denmin, den[[i]]$y)
        if (xx == "TML" && x[[i]]$call$errors == "logWeibull") 
            q[[i]] <- log(qweibull(ppoints(length(r[[i]])), shape = 1))
        else q[[i]] <- qnorm(ppoints(length(r[[i]])))
        y[[i]] <- f[[i]] + r[[i]]
        yname[[i]] <- deparse(form[[i]][[2]])
        fname[[i]] <- paste("Fitted :", deparse(form[[i]][[3]]),
            collapse = " ")
        if (xplots) {
            xx <- model.matrix(x[[i]])
            indep[[i]] <- xx
            for (j in 1:lnames) {
                k <- match(unames[j], dimnames(xx)[[2]], nomatch = 0)
                if (k) {
                  xmax[j] <- max(xmax[j], xx[, k])
                  xmin[j] <- min(xmin[j], xx[, k])
                }
            }
        }
    }
    minabs <- function(x) {
        min(abs(x), na.rm = TRUE)
    }
    maxabs <- function(x) {
        max(abs(x), na.rm = TRUE)
    }
    choose <- function(i) {
        if (i == 1)
            16
        else i
    }
    fmax <- max(sapply(f, max, na.rm = TRUE))
    fmin <- min(sapply(f, min, na.rm = TRUE))
    ymax <- max(sapply(y, max, na.rm = TRUE))
    ymin <- min(sapply(y, min, na.rm = TRUE))
    rmin <- min(sapply(r, min, na.rm = TRUE))
    rmax <- max(sapply(r, max, na.rm = TRUE))
    rsmin <- min(sapply(rs, min, na.rm = TRUE))
    rsmax <- max(sapply(rs, max, na.rm = TRUE))
    qmin <- min(sapply(q, min, na.rm = TRUE))
    qmax <- max(sapply(q, max, na.rm = TRUE))
    rrmin <- rmin - 0.07 * (rmax - rmin)
    rrmax <- rmax + 0.07 * (rmax - rmin)

    if(length(which)==1 & xplots==FALSE) ask <- FALSE
    if(show[1]){
        zz <- hist(unlist(list(r[[1]], c(rrmin, rrmax))), plot = FALSE,
            nclass = 10)
        breaks <- zz$breaks
        par(mfcol = c(1, length(x)))
        for (i in 1:length(x)) {
            hist(r[[i]], probability = TRUE, xlab = "Residuals", ylab = "Density",
                main = names[i], ylim = c(denmin, denmax), density = -1,
                breaks = breaks, cex.main = 1, cex.lab = 1, cex.axis = 1)
        }
    }
    par(xpd = FALSE)
    par(mfcol = c(1, 1))
    if(show[2]){
        plot(q[[1]], rs[[1]],type="n", xlab = "Theoretical Quantiles", ylab = "Standardized Residuals", main = "", xlim = c(qmin,
                qmax), ylim = c(rsmin, rsmax))
        for (i in 1:length(x)) {
            if (i != 1)
                par(new = TRUE)
            points(q[[i]], sort(rs[[i]]), pch = choose(i))
        }
        abline(0,1)
        legend(leg.position[1] , legend = names, cex = 0.8, bty = "n",
            pch = c(16, 2:length(x)))
    }
    if(show[3]){
        plot(f[[1]], y[[1]], xlab = "Fitted Values", ylab = "Response",
            pch = 16, xlim = c(fmin, fmax), ylim = c(ymin, ymax))
        abline(0, 1, lty = 2)
        for (i in 2:length(x)) {
            points(f[[i]], y[[i]], pch = i)
        }
        legend(leg.position[2], legend = names, bty = "n", cex = 0.8, pch = c(16, 2:length(x)))
    }
    if(show[4]){
        plot(f[[1]], r[[1]], xlab = "Fitted Values", ylab = "Residuals",
            pch = 16, xlim = c(fmin, fmax), ylim = c(rmin, rmax))
        for (i in 2:length(x)) {
            points(f[[i]], r[[i]], pch = i)
        }
        legend(leg.position[3], legend = names, bty = "n", cex = 0.8, pch = c(16, 2:length(x)))
    }
    if (xplots) {
        for (j in 1:length(unames)) {
            plt <- TRUE
            xin <- 1:length(x)
            par(mar = c(5, 4, 4, 10) + 0.1, xpd = TRUE)
            for (i in 1:length(x)) {
                k <- match(unames[j], dimnames(indep[[i]])[[2]],
                  nomatch = 0)
                if (k) {
                    if (plt) {
                        plt <- FALSE
                        plot(indep[[i]][, k], r[[i]], xlab = unames[j],
                            ylab = "Residuals", pch = choose(i), xlim = c(xmin[j],
                            xmax[j]), ylim = c(rmin, rmax))
                    }
                    else points(indep[[i]][, k], r[[i]], pch = choose(i))
                }
                else xin[i] <- -xin[i]
            }
            legend(xmax[j] + 0.1, rmax, legend = names[xin], bty = "n",
                pch = c(16, 2:length(xin)))
        }
    }
    invisible()
}
