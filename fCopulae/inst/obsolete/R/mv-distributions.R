
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

# Copyrights (C)
# for this R-port:
#   1999 - 2007, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:           PARAMETER ESTIMATION:
#  fMV                 S4 Object of class 'fMV'
#  mvFit               Fits a MV Normal or Student-t Distribution
#  print.fMV           S3: Print method for objects of class 'fMV'
#  plot.fMV            S3: Plot method for objects of class 'fMV'
#  summary.fMV         S3: Summary method for objects of class 'fMV'
#  .mvnormFit         Fits a Multivariate Normal Distribution
#  .mvstFit            Fits a Multivariate Student-t Distribution
#  .mvsnormPlot        Plots for Multivariate Normal Distribution
#  .mvstPlot           Plots for Multivariate Student-t Distribution
# REQUIREMENTS:       DESCRIPTION:
#  "mvtnorm"           Contributed R - Package
#  "sn" | "mnormt"     Contributed R - Package
################################################################################


################################################################################
# PARAMETER FIT:


setClass("fMV",
    representation(
        call = "call",
        method = "character",
        model = "list",
        data = "data.frame",
        fit = "list",
        title = "character",
        description = "character")
)


# ------------------------------------------------------------------------------


mvFit =
function(x, method = c("snorm", "st"), fixed.df = NA, title = NULL,
description = NULL, trace = FALSE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:

    # FUNCTION:

    # Fit:
    if (method[1] == "snorm") {
        # Normal Fit:
        fit = .mvsnormFit(x = x, trace = trace, ...)
        fit$df = Inf
    }
    if (method[1] == "st") {
       # Student-t Fit:
       fit = .mvstFit(x = x, fixed.df = fixed.df, trace = trace, ...)
    }

    # Add to fit:
    fit$method = method[1]
    class(fit) = "list"

    # Model Slot:
    model = list(beta = fit$beta, Omega = fit$Omega,
        alpha = fit$alpha, df = fit$df)

    # Title Slot:
    if (is.null(title)) {
        if (method[1] == "snorm")
            title = "Multivariate Normal Distribution"
        if (method[1] == "st")
            title = "Multivariate Student-t Distribution"
    }

    # Description Slot:
    if (is.null(description)) description = description()

    # Return Value:
    new("fMV",
        call = as.call(match.call()),
        method = as.character(method[1]),
        model = model,
        data = as.data.frame(x),
        fit = fit,
        title = as.character(title),
        description = as.character(description) )
}


# ------------------------------------------------------------------------------



setMethod("show", "fMV",
    function(object)
{   # A function implemented by Diethelm Wuertz

    # Description:

    # Arguments:

    # FUNCTION:

    # Extract fit:
    fit = object@fit

    # Print:
    cat("\nCall:\n ")
    print.default(fit$call)

    cat("\nParameter Sstimates:\n")
    print.default(fit$dp)

    cat("\nParameter Errors:\n")
    print.default(fit$se)

    # cat("\nOptimization:\n")
    # print.default(fit$optim)
})


# ------------------------------------------------------------------------------


plot.fMV =
function(x, which = "ask", ...)
{   # A function implemented by Diethelm Wuertz

    # Description:

    # Arguments:

    # FUNCTION:

    # Plot:
    if (x@fit$method == "snorm") {
        # Multivariate Skew Normal Distribution:
        return(.mvsnormPlot(x = x@fit, which = which, ...))
    }
    if (x@fit$method == "st") {
         # Multivariate Skew Student-t Distribution:
         return(.mvstPlot(x = x@fit, which = which, ...))
    }
}


# ------------------------------------------------------------------------------


summary.fMV =
function(object, which = "ask", doplot = TRUE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:

    # Arguments:

    # FUNCTION:

    # Print:
    print(x = object, ...)

    # Plot:
    if (doplot) plot(x = object, which = which, doplot, ...)

    # Return Value:
    invisible(object)
}


################################################################################
# INERNAL FUNCTIONS:


.mvsnormFit =
function(x, trace = FALSE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal Function

    # Arguments:

    # FUNCTION:

    # Settings:
    y = x
    y.name = deparse(substitute(y))
    y.names = dimnames(y)[[2]]
    y = as.matrix(y)
    colnames(y) = y.names
    k = ncol(y)
    freq = rep(1, nrow(y))
    n = sum(freq)
    X = rep(1, nrow(y))
    X = as.matrix(X)
    m = ncol(X)
    dimnames(y) = list(NULL, outer("V", as.character(1:k), paste, sep = ""))
    y.names = as.vector(dimnames(y)[[2]])
    qrX = qr(X)

    # Fit:
    mle = msn.mle(X = X, y = y, freq = freq, trace = trace, ...)
    mle$call = match.call()
    mle$y = y
    mle$y.names = y.names

    # Parameters:
    mle$beta = beta = mle$dp$beta
    mle$xi = xi = X %*% beta
    mle$Omega = Omega = mle$dp$Omega
    mle$alpha = alpha = mle$dp$alpha

    # Test:
    # dev.norm = msn.dev(c(qr.coef(qrX, y), rep(0, k)), X, y, freq)
    # test = dev.norm + 2 * mle$logL
    # p.value = 1 - pchisq(test, k)
    # mle$test.normality = list(LRT = test, p.value = p.value)

    # Save for Plot:
    Xb = qr.fitted(qrX, y)
    res = qr.resid(qrX, y)
    mle$k = k
    mle$n = n
    mle$pp = qchisq((1:n)/(n + 1), k)
    mle$rad.n = apply((y - Xb) * ((y - Xb) %*% solve(var(res))), 1, sum)
    mle$rad.sn = apply((y - xi) * ((y - xi) %*% solve(Omega)), 1, sum)

    # Return Value:
    class(mle) = "snFit"
    mle
}


# ------------------------------------------------------------------------------


.mvstFit =
function(x, fixed.df = NA, trace = FALSE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal Function

    # Arguments:

    # FUNCTION:

    # Settings:
    y = as.matrix(x)
    k = ncol(y)
    y.name = deparse(substitute(y))
    dimnames(y) = list(NULL, outer("V", as.character(1:k), paste, sep = ""))
    y.names = dimnames(y)[[2]]

    freq = rep(1, nrow(y))
    n = sum(freq)

    X = as.matrix(rep(1, nrow(y)))
    qrX = qr(X)
    m = ncol(X)

    # Fit:
    mle = mst.mle(X = X, y = y, freq = freq, fixed.df = fixed.df,
        trace = trace, ...)
    mle$call = match.call()
    mle$y = y
    mle$y.names = y.names

    # Parameters:
    mle$beta = beta = mle$dp$beta
    mle$xi = xi = X %*% beta
    mle$Omega = Omega = mle$dp$Omega
    mle$alpha = alpha = mle$dp$alpha
    mle$df = df = mle$dp$df

    # Save for Plot:
    Xb = qr.fitted(qrX, y)
    res = qr.resid(qrX, y)
    mle$k = k
    mle$n = n
    mle$pp = k * qf((1:n)/(n + 1), k, df)
    mle$rad.n = as.vector(apply(res * (res %*% solve(var(res))), 1, sum))
    mle$rad.sn = as.vector(apply((y - xi)*((y - xi) %*% solve(Omega)), 1, sum))

    # Return Value:
    class(mle) = "stFit"
    mle
}


# ------------------------------------------------------------------------------


.mvsnormPlot =
function(x, which = "ask", ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal Plot Function

    # Arguments:
    #   x - the slot @fit from an object of class "fMV"

    # FUNCTION:

    # Settings:
    dim = ncol(x$y)

    # Plot Title:
    plot1Title = "Scatterplots"
    if (dim == 1) plot1Title = "Histogram Plot"

    # Plot:
    interactivePlot(
        x = x,
        choices = c(
            plot1Title,
            "Normal QQ-Plot",
            "Skew-Normal QQ-Plot",
            "Normal PP-Plot",
            "Skew-Normal PP-Plot"),
        plotFUN = c(
            ".mvsnorm.plot.1",
            ".mvsnorm.plot.2",
            ".mvsnorm.plot.3",
            ".mvsnorm.plot.4",
            ".mvsnorm.plot.5"),
        which = which)

    # Return Value:
    invisible(x)
}


# ------------------------------------------------------------------------------


.mvsnorm.plot.1 <-
function(x)
{
    # Plot:
    dim = x$k
    if(dim == 1) .mvsnorm.plot.1A(x) else .mvsnorm.plot.1B(x)
}


# ------------------------------------------------------------------------------


.mvsnorm.plot.1A <-
function(x)
{
    # Plot:
    z = x
    y0 <- z$y
    xi0 <- apply(z$xi, 2, mean)
    y0 <- as.vector(y0)
    x <- seq(min(pretty(y0, 10)), max(pretty(y0, 10)), length = 100)
    omega <- sqrt(diag(z$Omega))
    dp0 <- c(xi0, omega, z$alpha)
    xlab <- z$y.name
    hist(y0, prob = TRUE, breaks = "FD", xlab = xlab,
        ylab = "density", border = "white", col = "steelblue4",
        main = z$y.name)
    lines(x, dsn(x, dp0[1], dp0[2], dp0[3]))
    if (length(y0) < 201)
        points(y0, rep(0, z$n), pch = 1)
}


# ------------------------------------------------------------------------------


.mvsnorm.plot.1B <-
function(x)
{
    # Plot:
    opt = options()
    options(warn = -1)
    pairs(
        x$y,
        labels = x$y.names,
        panel = function(x, y, Y, y.names, xi, Omega, alpha) {
            for (i in 1:length(alpha)) {
                if (all(Y[, i] == x))
                    Ix = i
                if (all(Y[, i] == y))
                    Iy = i }
            points(x, y)
            marg = msn.marginal(xi, Omega, alpha, c(Ix, Iy))
            xi.marg = marg$xi
            Omega.marg = marg$Omega
            alpha.marg = marg$alpha
            x1 = seq(min(x), max(x), length = 30)
            x2 = seq(min(y), max(y), length = 30)
            dsn2.plot(x1, x2, xi.marg, Omega.marg, alpha.marg,
                add = TRUE, col = "steelblue4")},
        Y = x$y,
        y.names = dimnames(x$y)[[2]],
        xi = apply(x$xi, 2, mean),
        Omega = x$Omega,
        alpha = x$alpha)
    options(opt)
}


# ------------------------------------------------------------------------------


.mvsnorm.plot.2 <-
function(x)
{
    # Plot:
    plot(x$pp, sort(x$rad.n), pch = 1, ylim = c(0, max(x$rad.n, x$rad.sn)),
        xlab = "Chi-square Percentiles",
        ylab = "Mahalanobis Distances")
    abline(0, 1, lty = 3)
    title(main = "Normal QQ-Plot", sub = x$y.name)
}


# ------------------------------------------------------------------------------


.mvsnorm.plot.3 <-
function(x)
{
    # Plot:
    plot(x$pp, sort(x$rad.sn), pch = 1, ylim = c(0, max(x$rad.n, x$rad.sn)),
        xlab = "Percentiles of chi-square distribution",
        ylab = "Mahalanobis distances")
    abline(0, 1, lty = 3)
    title(main = "Skew-Normal QQ-Plot", sub = x$y.name)
}


# ------------------------------------------------------------------------------


.mvsnorm.plot.4 <-
function(x)
{
    # Plot:
    plot((1:x$n)/(x$n + 1), sort(pchisq(x$rad.n, x$k)),
        xlab = "",  ylab = "")
    abline(0, 1, lty = 3)
    title(main = "Normal PP-Plot", sub = x$y.name)
}


# ------------------------------------------------------------------------------


.mvsnorm.plot.5 <-
function(x)
{
    # Plot:
    plot((1:x$n)/(x$n + 1), sort(pchisq(x$rad.sn, x$k)),
        xlab = "", ylab = "")
    abline(0, 1, lty = 3)
    title(main = "Skew-Normal PP-Plot", sub = x$y.name)
}


# ------------------------------------------------------------------------------


.mvstPlot =
function(x, which = "ask", ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal Plot Function

    # Arguments:
    #   x - the slot @fit from an object of class "fMV"

    # FUNCTION:

    # Settings:
    dim = ncol(x$y)

    # Plot Title:
    plot1Title = "Scatterplots"
    if (dim == 1) plot1Title = "Histogram Plot"

    # Plot:
    plot1Title = "Scatterplots"
    if (dim == 1) plot1Title = "Histogram Plot"
    interactivePlot(
        x = x,
        choices = c(
            plot1Title,
            "Normal QQ-Plot",
            "Skew-Normal QQ-Plot",
            "Normal PP-Plot",
            "Skew-Normal PP-Plot"),
        plotFUN = c(
            ".mvst.plot.1",
            ".mvst.plot.2",
            ".mvst.plot.3",
            ".mvst.plot.4",
            ".mvst.plot.5"),
        which = which)

    # Return Value:
    invisible(x)
}


# ------------------------------------------------------------------------------


.mvst.plot.1 <-
function(x)
{
    # Plot:
    dim = x$k
    if(dim == 1) .mvst.plot.1A(x) else .mvst.plot.1B(x)
}


# ------------------------------------------------------------------------------


.mvst.plot.1A <-
function(x)
{
    # Plot:
    z = x
    y0 <- z$y
    xi0 <- apply(z$xi, 2, mean)
    y0 <- as.vector(y0)
    x <- seq(min(pretty(y0, 10)), max(pretty(y0, 10)), length = 100)
    omega <- sqrt(diag(z$Omega))
    dp0 <- c(xi0, omega, z$alpha, z$df)
    xlab <- z$y.name
    hist(y0, prob = TRUE, breaks = "FD", xlab = xlab,
        ylab = "density", border = "white", col = "steelblue4",
        main = z$y.name)
    lines(x, dst(x, dp0[1], dp0[2], dp0[3], dp0[4]))
    if (length(y0) < 201)
        points(y0, rep(0, z$n), pch = 1)
}


# ------------------------------------------------------------------------------


.mvst.plot.1B <-
function(x)
{
    # Plot:
    opt = options()
    options(warn = -1)
    pairs(
        x$y,
        labels = x$y.names,
        panel = function(x, y, Y, y.names, xi, Omega, alpha, df) {
            for (i in 1:length(alpha)) {
                if (all(Y[, i] == x))
                    Ix = i
                if (all(Y[, i] == y))
                    Iy = i }
            points(x, y)
            marg = msn.marginal(xi, Omega, alpha, c(Ix, Iy))
            xi.marg = marg$xi
            Omega.marg = marg$Omega
            alpha.marg = marg$alpha
            x1 = seq(min(x), max(x), length = 30)
            x2 = seq(min(y), max(y), length = 30)
            dst2.plot(x1, x2, xi.marg, Omega.marg, alpha.marg,
                df, add = TRUE, col = "steelblue4")} ,
        Y = x$y,
        y.names = dimnames(x$y)[[2]],
        xi = apply(x$xi, 2, mean),
        Omega = x$Omega,
        alpha = x$alpha,
        df = x$df)
    options(opt)
}


# ------------------------------------------------------------------------------


.mvst.plot.2 <-
function(x)
{
    # Plot:
    plot(x$pp, sort(x$rad.n), pch = 1, ylim = c(0, max(x$rad.n, x$rad.sn)),
        xlab = "Chi-square Percentiles",
        ylab = "Mahalanobis Distances")
    abline(0, 1, lty = 3)
    title(main = "Normal QQ-Plot", sub = x$y.name)
}


# ------------------------------------------------------------------------------


.mvst.plot.3 <-
function(x)
{
    # Plot:
    plot(x$pp, sort(x$rad.sn), pch = 1, ylim = c(0, max(x$rad.n, x$rad.sn)),
        xlab = "Percentiles of chi-square distribution",
        ylab = "Mahalanobis distances")
    abline(0, 1, lty = 3)
    title(main = "Skew-Normal QQ-Plot", sub = x$y.name)
}


# ------------------------------------------------------------------------------


.mvst.plot.4 <-
function(x)
{
    # Plot:
    plot((1:x$n)/(x$n + 1), sort(pchisq(x$rad.n, x$k)),
        xlab = "",  ylab = "")
    abline(0, 1, lty = 3)
    title(main = "Normal PP-Plot", sub = x$y.name)
}


# ------------------------------------------------------------------------------


.mvst.plot.5 <-
function(x)
{
    # Plot:
    plot((1:x$n)/(x$n + 1), sort(pchisq(x$rad.sn, x$k)),
        xlab = "", ylab = "")
    abline(0, 1, lty = 3)
    title(main = "Skew-Normal PP-Plot", sub = x$y.name)
}


################################################################################

