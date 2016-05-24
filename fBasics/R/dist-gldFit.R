
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:            DESCRIPTION:
#  gldFit               Fits parameters of a GLD
#  .gldFit.mle          Fits parameters of a GLD using maximum log-likelihood
#  .gldFit.mps          Fits parameters of a GLD using maximum product spacings
#  .gldFit.gof          Fits parameters of a GLD using GoF Statistics
#   .ksGLD               Kolmogorov Smirnov Statistics
#   .cvmGLD              Cramer von Mise Statistics
#   .adGLD               Anderson Darling Statistics
#  .gldFit.hist         Fits parameters of a GLD using a histogram fit
#   type="fd"            Freedman-Diaconis binning
#   type="scott"         Scott binning
#   type="sturges"       Sturges binning
#  .gldFit.rob          Fits parameters of a GLD using robust moments fit
################################################################################


gldFit <-
function(x, lambda1 = 0, lambda2 = -1, lambda3 = -1/8, lambda4 = -1/8,
    method = c("mle", "mps", "gof", "hist", "rob"),
    scale = NA, doplot = TRUE, add = FALSE, span = "auto", trace = TRUE,
    title = NULL, description = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits parameters of a GLD using maximum log-likelihood

    # Arguments:
    #   lambda1 - a numeric value, the location parameter
    #   lambda2 - a numeric value, the scale parameter
    #   lambda3 - a numeric value, the first shape parameter
    #   lambda4 - a numeric value, the second parameter

    # Note:
    #   The GLD uses RS parameterization and parameters in Region 4,
    #   i.e. lambda[3,4]<0, lambda2<0

    # FUNCTION:

    # Check Parameters:
    stopifnot(lambda2 < 0)
    stopifnot(lambda3 < 0)
    stopifnot(lambda4 < 0)

    # Settings:
    method = match.arg(method)

    # Parameter Fit:
    if (method == "mle") {
        ans = .gldFit.mle(x, lambda1, lambda2, lambda3, lambda4,
            scale, doplot, add, span, trace, title, description, ...)
    } else if (method == "mps") {
        ans = .gldFit.mps(x, lambda1, lambda2, lambda3, lambda4,
            scale, doplot, add, span, trace, title, description, ...)
    } else if (method == "gof") {
        ans = .gldFit.gof(x, lambda1, lambda2, lambda3, lambda4,
            scale, doplot, add, span, trace, title, description, ...)
    } else if (method == "hist") {
        ans = .gldFit.hist(x, lambda1, lambda2, lambda3, lambda4,
            scale, doplot, add, span, trace, title, description, ...)
    } else if (method == "rob") {
        ans = .gldFit.rob (x, lambda1, lambda2, lambda3, lambda4,
            scale, doplot, add, span, trace, title, description, ...)
    }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.gldFit.mle <-
function(x, lambda1 = 0, lambda2 = -1, lambda3 = -1/8, lambda4 = -1/8,
    scale = NA, doplot = TRUE, add = FALSE, span = "auto", trace = TRUE,
    title = NULL, description = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits parameters of a GLD using maximum log-likelihood

    # Example:
    #   require(fBasics)
    #   set.seed(4711); x=rgld(5000); fit=gldFit.mle(x)@fit$estimate; fit

    # FUNCTION:

    # Settings:
    scale = FALSE

    # Transform:
    x.orig = x
    x = as.vector(x)
    if (scale) {
        SD = sd(x)
        x = x / SD }

    # Settings:
    CALL = match.call()

    # Objective Function:
    obj <- function(x, y = x, trace) {
        DGLD = try(dgld(y, x[1], x[2], x[3], x[4]), silent = TRUE)
        if (class(DGLD) == "try-error") return(1e9)
        f = -sum(log(DGLD))
        # Print Iteration Path:
        if (trace) {
            cat("\n Objective Function Value:  ", f)
            cat("\n Parameter Estimates:       ", x[1], x[2], x[3], x[4], "\n")
        }
        f }

    # Parameter Estimation:
    eps = 1e-10
    BIG = 100
    r = nlminb(
        start = c(lambda1, lambda2, lambda3, lambda4),
        objective = obj,
        lower = c(-BIG, -BIG, -BIG, -BIG),
        upper = c(+BIG, -eps, -eps, -eps),
        y = x,
        trace = trace)
    names(r$par) <- c("lambda1", "lambda2", "lambda3", "lambda4")

    # Add Title and Description:
    if (is.null(title)) title = "GLD Region 4 Parameter Estimation"
    if (is.null(description)) description = description()

    # Rescale Result:
    if (scale) {
        r$par = r$par / c(SD, SD, 1/SD, 1/SD)
        r$objective = obj(r$par, y = as.vector(x.orig), trace = trace)
    }
    fit = list(estimate = r$par, minimum = -r$objective, code = r$convergence)

    # Optional Plot:
    if (doplot) {
        x = as.vector(x.orig)
        if (span == "auto") span = seq(min(x), max(x), length = 501)
        z = density(x, ...)
        x = z$x[z$y > 0]
        y = z$y[z$y > 0]
        y.points = dgld(span, r$par[1], r$par[2], r$par[3], r$par[4])
        ylim = log(c(min(y.points), max(y.points)))
        if (add) {
            lines(x = span, y = log(y.points), ...)
        } else {
            plot(x, log(y), xlim = c(span[1], span[length(span)]),
                ylim = ylim, type = "p", xlab = "x", ylab = "log f(x)", ...)
            title(main = title)
            lines(x = span, y = log(y.points), col = "steelblue")
        }
    }

    # Return Value:
    new("fDISTFIT",
        call = as.call(CALL),
        model = "GLD Region 4 Distribution",
        data = as.data.frame(x.orig),
        fit = fit,
        title = as.character(title),
        description = description() )
}


# ------------------------------------------------------------------------------


.gldFit.mps <-
function(x, lambda1 = 0, lambda2 = -1, lambda3 = -1/8, lambda4 = -1/8,
    type = c("sum", "mean", "max", "median", "var"),
    scale = NA, doplot = TRUE, add = FALSE, span = "auto", trace = TRUE,
    title = NULL, description = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits parameters of a GLD using maximum product spacings

    # Example:
    #   require(fBasics)
    #   set.seed(4711); x=rgld(5000); fit=.gldFit.mps(x)@fit$estimate; fit

    # FUNCTION:

    # Settings:
    scale = FALSE
    type = match.arg(type)
    CALL = match.call()

    # Transform:
    x.orig = x
    x = as.vector(x)
    if (scale) {
        SD = sd(x)
        x = x / SD }

    # Objective Function:
    TYPE = toupper(type)
    if(type == "sum") {
        typeFun = sum
    } else if (type == "mean") {
        typeFun = mean
    } else if (type == "median") {
        typeFun = median
    } else if (type == "max") {
        typeFun = function(x) -max(x)
    } else if (type == "var") {
        typeFun = function(x) -var(x) }
    obj = function(x, y = x, typeFun, trace) {
        PGLD = try(pgld(sort(y), x[1], x[2], x[3], x[4]), silent = TRUE)
        if (class(PGLD) == "try-error") return(1e9)
        DH = diff(c(0, na.omit(PGLD), 1))
        f = try(-typeFun(log(DH[DH > 0])), silent = TRUE)
        if (class(PGLD) == "try-error") return(1e9)
        # Print Iteration Path:
        if (trace) {
            cat("\n Objective Function Value:  ", f)
            cat("\n Parameter Estimates:       ", x[1], x[2], x[3], x[4], "\n")
        }
        f }

    # Parameter Estimation:
    eps = 1e-10
    BIG = 100
    r = nlminb(
        start = c(lambda1, lambda2, lambda3, lambda4),
        objective = obj,
        lower = c(-BIG, -BIG, -BIG, -BIG),
        upper = c(+BIG, -eps, -eps, -eps),
        y = x,
        typeFun = typeFun,
        trace = trace)
    names(r$par) <- c("lambda1", "lambda2", "lambda3", "lambda4")

    # Add Title and Description:
    if (is.null(title)) title = "GLD Region 4 MPS Estimation"
    if (is.null(description)) description = description()

    # Rescale Result:
    if (scale) {
        r$par = r$par / c(SD, SD, 1/SD, 1/SD)
        r$objective = obj(r$par, y = as.vector(x.orig), trace = trace)
    }
    fit = list(estimate = r$par, minimum = -r$objective, code = r$convergence)

    # Optional Plot:
    if (doplot) {
        x = as.vector(x.orig)
        if (span == "auto") span = seq(min(x), max(x), length = 501)
        z = density(x, ...)
        x = z$x[z$y > 0]
        y = z$y[z$y > 0]
        y.points = dgld(span, r$par[1], r$par[2], r$par[3], r$par[4])
        ylim = log(c(min(y.points), max(y.points)))
        if (add) {
            lines(x = span, y = log(y.points), ...)
        } else {
            plot(x, log(y), xlim = c(span[1], span[length(span)]),
                ylim = ylim, type = "p", xlab = "x", ylab = "log f(x)", ...)
            title(main = title)
            lines(x = span, y = log(y.points), col = "steelblue")
        }
    }

    # Return Value:
    new("fDISTFIT",
        call = as.call(CALL),
        model = "GLD Region 4 Distribution",
        data = as.data.frame(x.orig),
        fit = fit,
        title = as.character(title),
        description = description() )
}


# ------------------------------------------------------------------------------


.ksGLD <-
function(N, PFGL) {
    D = 1/(2*N)  + max ( abs( PFGL - ((1:N)-0.5)/N ) )
    D }

.cvmGLD <-
function(N, PFGL) {
    W2 = 1/(12*N)  + sum ( PFGL - ((1:N)-0.5)/N  )^2
    W2 }

.adGLD <-
function(N, PFGL) {
    A2 = -N -(1/N) * sum( (2*(1:N)-1) * ( log(PFGL) + log(1-rev(PFGL)) ) )
    A2 }


# ------------------------------------------------------------------------------


.gldFit.gof <-
function(x, lambda1 = 0, lambda2 = -1, lambda3 = -1/8, lambda4 = -1/8,
    type = c("ad", "ks", "cvm"),
    scale = NA, doplot = TRUE, add = FALSE, span = "auto", trace = TRUE,
    title = NULL, description = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits parameters of a GLD using GoF Statistics

    # Example:
    #   require(fBasics)
    #   set.seed(4711); x=rgld(5000); fit=.gldFit.gof(x)@fit$estimate; fit

    # FUNCTION:

    # Settings:
    scale = FALSE
    type = match.arg(type)
    CALL = match.call()

    # Transform:
    x.orig = x
    x = as.vector(x)
    if (scale) {
        SD = sd(x)
        x = x / SD }

    # Objective Function:
    TYPE = toupper(type)
    ## DJS 20/02/2010
    ## typeFun <- match.fun(paste(".", type, sep = ""))
    typeFun <- match.fun(paste(".", type, "GLD", sep = ""))
    obj = function(x, y = x, typeFun, trace) {
        PFGL = try(pgld(sort(y), x[1], x[2], x[3], x[4]), silent = TRUE)
        if (class(PFGL) == "try-error") return(1e9)
        PFGL = PFGL[PFGL > 0]
        PFGL = PFGL[PFGL < 1]
        N = length(PFGL)
        f = typeFun(N, PFGL)
        if (is.na(f)) return(1e9)
        # Print Iteration Path:
        if (trace) {
            cat("\n Objective Function Value:  ", f)
            cat("\n Parameter Estimates:       ", x[1], x[2], x[3], x[4],"\n")
        }
        f }

    # Parameter Estimation:
    eps = 1e-10
    BIG = 100
    r = nlminb(
        start = c(lambda1, lambda2, lambda3, lambda4),
        objective = obj,
        lower = c(-BIG, -BIG, -BIG, -BIG),
        upper = c(+BIG, -eps, -eps, -eps),
        y = x,
        typeFun = typeFun,
        trace = trace)
    names(r$par) <- c("lambda1", "lambda2", "lambda3", "lambda4")

    # Add Title and Description:
    if (is.null(title)) title = paste("GLD Region 4", TYPE, "Estimation")
    if (is.null(description)) description = description()

    # Rescale Result:
    if (scale) {
        r$par = r$par / c(SD, SD, 1/SD, 1/SD)
        r$objective = obj(r$par, y = as.vector(x.orig), trace = trace)
    }
    fit = list(estimate = r$par, minimum = -r$objective, code = r$convergence)

    # Optional Plot:
    if (doplot) {
        x = as.vector(x.orig)
        if (span == "auto") span = seq(min(x), max(x), length = 501)
        z = density(x, ...)
        x = z$x[z$y > 0]
        y = z$y[z$y > 0]
        y.points = dgld(span, r$par[1], r$par[2], r$par[3], r$par[4])
        ylim = log(c(min(y.points), max(y.points)))
        if (add) {
            lines(x = span, y = log(y.points), ...)
        } else {
            plot(x, log(y), xlim = c(span[1], span[length(span)]),
                ylim = ylim, type = "p", xlab = "x", ylab = "log f(x)", ...)
            title(main = title)
            lines(x = span, y = log(y.points), col = "steelblue")
        }
    }

    # Return Value:
    new("fDISTFIT",
        call = as.call(CALL),
        model = "GLD Region 4 Distribution",
        data = as.data.frame(x.orig),
        fit = fit,
        title = as.character(title),
        description = description() )
}


# ------------------------------------------------------------------------------


.gldFit.hist <-
function(x, lambda1 = 0, lambda2 = -1, lambda3 = -1/8, lambda4 = -1/8,
    type = c("fd", "sturges", "scott"),
    scale = NA, doplot = TRUE, add = FALSE, span = "auto", trace = TRUE,
    title = NULL, description = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits parameters of a GLD using a histogram fit

    # Example:
    #   require(fBasics)
    #   set.seed(4711); x=rgld(5000); fit=gldFit.hist(x)@fit$estimate; fit

    # FUNCTION:

    # Settings:
    scale = FALSE
    CALL = match.call()
    type = match.arg(type)
    if (type == "fd") type = "FD"
    else if (type == "scott") type = "Scott"
    else if (type == "sturges") type = "Sturges"

    # Transform:
    x.orig = x
    x = as.vector(x)
    if (scale) {
        SD = sd(x)
        x = x / SD }

    # Histogram:
    HIST = hist(x, breaks = type, plot = FALSE)

    # Objective Function:
    obj = function(x, hist, trace) {
        DFGL = try(dgld(hist$mids, x[1], x[2], x[3], x[4]), silent = TRUE)
        if (class(DFGL) == "try-error") return(1e9)
        DST = (hist$density - DFGL)^2
        f = sum(DST)
        # Print Iteration Path:
        if (trace) {
            cat("\n Objective Function Value:  ", f)
            cat("\n Parameter Estimates:       ", x[1], x[2], x[3], x[4], "\n")
        }
        f }

    # Parameter Estimation:
    eps = 1e-10
    BIG = 100
    r = nlminb(
        start = c(lambda1, lambda2, lambda3, lambda4),
        objective = obj,
        lower = c(-BIG, -BIG, -BIG, -BIG),
        upper = c(+BIG, -eps, -eps, -eps),
        hist = HIST,
        trace = trace)
    names(r$par) <- c("lambda1", "lambda2", "lambda3", "lambda4")

    # Add Title and Description:
    if (is.null(title)) title = "GLD Region 4 Histogram Estimation"
    if (is.null(description)) description = description()

    # Rescale Result:
    if (scale) {
        r$par = r$par / c(SD, SD, 1/SD, 1/SD)
        r$objective = obj(r$par, y = as.vector(x.orig), trace = trace)
    }
    fit = list(estimate = r$par, minimum = -r$objective, code = r$convergence)

    # Optional Plot:
    if (doplot) {
        x = as.vector(x.orig)
        if (span == "auto") span = seq(min(x), max(x), length = 501)
        z = density(x, ...)
        x = z$x[z$y > 0]
        y = z$y[z$y > 0]
        y.points = dgld(span, r$par[1], r$par[2], r$par[3], r$par[4])
        ylim = log(c(min(y.points), max(y.points)))
        if (add) {
            lines(x = span, y = log(y.points), ...)
        } else {
            plot(x, log(y), xlim = c(span[1], span[length(span)]),
                ylim = ylim, type = "p", xlab = "x", ylab = "log f(x)", ...)
            title(main = title)
            lines(x = span, y = log(y.points), col = "steelblue")
        }
    }

    # Return Value:
    new("fDISTFIT",
        call = as.call(CALL),
        model = "GLD Region 4 Distribution",
        data = as.data.frame(x.orig),
        fit = fit,
        title = as.character(title),
        description = description() )
}


# ------------------------------------------------------------------------------


.gldFit.rob <-
function(x, lambda1 = 0, lambda2 = -1, lambda3 = -1/8, lambda4 = -1/8,
    scale = NA, doplot = TRUE, add = FALSE, span = "auto", trace = TRUE,
    title = NULL, description = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits parameters of a GLD using robust moments (quantile) fit

    # Example:
    #   require(fBasics)
    #   set.seed(4711); x=rgld(5000); fit=.gldFit.rob(x)@fit$estimate; fit

    # FUNCTION:

    # Settings:
    scale = FALSE

    # Transform:
    x.orig = x
    x = as.vector(x)
    if (scale) {
        SD = sd(x)
        x = x / SD }

    # Settings:
    CALL = match.call()
    xMED = sampleMED(x)
    xIQR = sampleIQR(x)
    xSKEW = sampleSKEW(x)
    xKURT = sampleKURT(x)

    # Objective Function:
    obj = function(x, xMED, xIQR, xSKEW, xKURT, trace)
    {
        lambda1 = x[1]
        lambda2 = x[2]
        lambda3 = x[3]
        lambda4 = x[4]
        ROB = try(sqrt(
            (xMED -gldMED (lambda1, lambda2, lambda3, lambda4))^2 +
            (xIQR -gldIQR (lambda1, lambda2, lambda3, lambda4))^2 +
            (xSKEW-gldSKEW(lambda1, lambda2, lambda3, lambda4))^2 +
            (xKURT-gldKURT(lambda1, lambda2, lambda3, lambda4))^2 ),
            silent = TRUE)
        if(class(ROB) == "try-error") return(1e9)
        f = ROB
        # Print Iteration Path:
        if (trace) {
            cat("\n Objective Function Value:  ", f)
            cat("\n Parameter Estimates:       ", x[1], x[2], x[3], x[4], "\n")
        }
        f
    }

    # Parameter Estimation:
    eps = 1e-10
    BIG = 100
    r = nlminb(
        start = c(lambda1, lambda2, lambda3, lambda4),
        objective = obj,
        lower = c(-BIG, -BIG, -BIG, -BIG),
        upper = c(+BIG, -eps, -eps, -eps),
        xMED = xMED, xIQR = xIQR, xSKEW = xSKEW, xKURT = xKURT,
        trace = trace)
    names(r$par) <- c("lambda1", "lambda2", "lambda3", "lambda4")

    # Add Title and Description:
    if (is.null(title)) title = "GLD Region 4 Robust Moment Estimation"
    if (is.null(description)) description = description()

    # Rescale Result:
    if (scale) {
        r$par = r$par / c(SD, SD, 1/SD, 1/SD)
        r$objective = obj(r$par, y = as.vector(x.orig), trace = trace)
    }
    fit = list(estimate = r$par, minimum = -r$objective, code = r$convergence)

    # Optional Plot:
    if (doplot) {
        x = as.vector(x.orig)
        if (span == "auto") span = seq(min(x), max(x), length = 501)
        z = density(x, ...)
        x = z$x[z$y > 0]
        y = z$y[z$y > 0]
        y.points = dgld(span, r$par[1], r$par[2], r$par[3], r$par[4])
        ylim = log(c(min(y.points), max(y.points)))
        if (add) {
            lines(x = span, y = log(y.points), ...)
        } else {
            plot(x, log(y), xlim = c(span[1], span[length(span)]),
                ylim = ylim, type = "p", xlab = "x", ylab = "log f(x)", ...)
            title(main = title)
            lines(x = span, y = log(y.points), col = "steelblue")
        }
    }

    # Return Value:
    new("fDISTFIT",
        call = as.call(CALL),
        model = "GLD Region 4 Distribution",
        data = as.data.frame(x.orig),
        fit = fit,
        title = as.character(title),
        description = description() )
}


################################################################################

