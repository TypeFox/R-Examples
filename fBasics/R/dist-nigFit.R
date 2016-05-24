
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
# FUNCTION:            GENERALIZED DISTRIBUTION:
#  nigFit               Fits parameters of a normal inverse Gaussian density
#  .nigFit.mle          max Log-likelihood Estimation
#  .nigFit.gmm          gemeralized method of moments estimation
#  .nigFit.mps          maximum product spacings estimation
#  .nigFit.vmps         minimum variance product spacings estimation
################################################################################


nigFit <- 
function(x, alpha = 1, beta = 0, delta = 1, mu = 0, 
    method = c("mle", "gmm", "mps", "vmps"), 
    scale = TRUE, doplot = TRUE, span = "auto", trace = TRUE, 
    title = NULL, description = NULL, ...)
{
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION: 
    
    # Settings:
    method = match.arg(method)
    
    # Select:
    if (method == "mle") {
        # MLE:
        fit = .nigFit.mle(x = x, alpha = alpha, beta = beta, delta = delta, 
            mu = mu , scale = scale, doplot = doplot, span = span, 
            trace = trace, title = title, description = description, ...)
    } else if (method == "gmm") {
        # GMM:
        fit = .nigFit.gmm(x = x, alpha = alpha, beta = beta, delta = delta, 
            mu = mu , scale = scale, doplot = doplot, span = span, 
            trace = trace, title = title, description = description, ...)
    } else if (method == "mps") {
        # MPS:
        fit = .nigFit.mps(x = x, alpha = alpha, beta = beta, delta = delta, 
            mu = mu , scale = scale, doplot = doplot, span = span, 
            trace = trace, title = title, description = description, ...)
    } else if (method == "vmps") {
        # MPS:
        fit = .nigFit.vmps(x = x, alpha = alpha, beta = beta, delta = delta, 
            mu = mu , scale = scale, doplot = doplot, span = span, 
            trace = trace, title = title, description = description, ...)
    }  
    
    # Return Value:
    fit
}


# ------------------------------------------------------------------------------


.nigFit.mle <-
function(x, alpha = 1, beta = 0, delta = 1, mu = 0, 
    scale = TRUE, doplot = TRUE, add = FALSE, span = "auto", trace = TRUE,
    title = NULL, description = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits parameters of a NIG using maximum log-likelihood  

    # Example:
    #   set.seed(4711); x = rnig(500); mle = .nigFit.mle(x); mle@fit$estimate
    
    # FUNCTION:

    # Transform:
    x.orig = x
    x = as.vector(x)
    if (scale) {
        SD = sd(x)
        x = x / SD }

    # Objective Function:
    obj = function(x, y = x, trace) {
        if (abs(x[2]) >= x[1]) return(1e9)
        f = -sum(dnig(y, x[1], x[2], x[3], x[4], log = TRUE))
        # Print Iteration Path:
        if (trace) {
            cat("\n Objective Function Value:  ", -f)
            cat("\n Parameter Estimates:       ", x, "\n")
        }
        f }
        
    # Parameter Estimation:
    eps = 1e-10
    BIG = 1000
    fit = nlminb(
        start = c(alpha, beta, delta, mu), 
        objective = obj,
        lower = c(eps, -BIG, eps, -BIG), 
        upper = BIG, 
        y = x, 
        trace = trace)
    names(fit$par) <- c("alpha", "beta", "delta", "mu")
    # Rescale Result:
    if (scale) {
        fit$scaleParams = c(SD, SD, 1/SD, 1/SD)
        fit$par = fit$par / fit$scaleParams
        fit$objective = obj(fit$par, y = as.vector(x.orig), trace = trace)
    } else {
        fit$scaleParams = rep(1, time = length(fit$par))
    }
    fit$scale = scale
    fit$estimate = fit$par
    fit$minimum = -fit$objective
    fit$code = fit$convergence
       
    # Standard Errors and t-Values:
    fit = .distStandardErrors(fit, obj, x)

    # Add Title and Description:
    if (is.null(title)) title = "Normal Inverse Gaussian Parameter Estimation"
    if (is.null(description)) description = description()

    # Optional Plot:
    if (doplot) .distFitPlot(
        fit, 
        x = x.orig, 
        FUN = "dnig", 
        main = "NIG Parameter Estimation", 
        span = span, add = add, ...)

    # Return Value:
    new("fDISTFIT",
        call = match.call(),
        model = "Normal Inverse Gaussian Distribution",
        data = as.data.frame(x.orig),
        fit = fit,
        title = as.character(title),
        description = description() )
}


# ------------------------------------------------------------------------------


.nigFit.gmm <-
function(x, 
    scale = TRUE, doplot = TRUE, add = FALSE, span = "auto", trace = TRUE,
    title = NULL, description = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits parameters of a NIG using GMM estimator  

    # Example:
    #   set.seed(4711); x = rnig(500); gmm = .nigFit.gmm(x)@fit$estimate; gmm
    
    # FUNCTION:

    # Transform:
    x.orig = x
    x = as.vector(x)
    if (scale) {
        SD = sd(x)
        x = x / SD }

    # Settings:
    CALL = match.call()

    # Parameter Estimation:
    obj <- function(Theta, x) {
        # Parameters:
        alpha = Theta[1]
        beta  = Theta[2]
        delta = Theta[3]
        mu    = Theta[4] 
        names(Theta) <- c("alpha", "beta", "delta", "mu")
        # Trace:
        if (TRUE) print(Theta)
        # Moments:
        m1 <- x   - .ghMuMoments(1, alpha, beta, delta, mu, lambda = -0.5)
        m2 <- x^2 - .ghMuMoments(2, alpha, beta, delta, mu, lambda = -0.5)
        m3 <- x^3 - .ghMuMoments(3, alpha, beta, delta, mu, lambda = -0.5)
        m4 <- x^4 - .ghMuMoments(4, alpha, beta, delta, mu, lambda = -0.5)
        # Result:
        f <- cbind(m1, m2, m3, m4)
        return(f)
    }
    r <- .gmm(g = obj, x = x, t0 = c(1, 0, 1, 0)) 
    names(r$par) <- c("alpha", "beta", "delta", "mu")
    
    # Add Title and Description:
    if (is.null(title)) title = "Normal Inverse Gaussian Parameter Estimation"
    if (is.null(description)) description = description()

    # Rescale Result:
    if (scale) {
        r$par = r$par / c(SD, SD, 1/SD, 1/SD)
        r$objective = NA
    }
    fit = list(estimate = r$par)

    # Optional Plot:
    if (doplot) {
        x = as.vector(x.orig)
        if (span == "auto") span = seq(min(x), max(x), length = 51)
        z = density(x, n = 100, ...)
        x = z$x[z$y > 0]
        y = z$y[z$y > 0]
        y.points = dnig(span, r$par[1], r$par[2], r$par[3], r$par[4])
        ylim = log(c(min(y.points), max(y.points)))
        if (add) {
            lines(x = span, y = log(y.points), col = "steelblue")
        } else {
            plot(x, log(y), xlim = c(span[1], span[length(span)]),
                ylim = ylim, type = "p", xlab = "x", ylab = "log f(x)", ...)
            title("NIG GMM Parameter Estimation")
            lines(x = span, y = log(y.points), col = "steelblue")
        }
    }

    # Return Value:
    new("fDISTFIT",
        call = as.call(CALL),
        model = "Normal Inverse Gaussian Distribution",
        data = as.data.frame(x.orig),
        fit = fit,
        title = as.character(title),
        description = description() )
}


# ------------------------------------------------------------------------------


.nigFit.mps <-
function(x, alpha = 1, beta = 0, delta = 1, mu = 0, 
    scale = TRUE, doplot = TRUE, add = FALSE, span = "auto", trace = TRUE,
    title = NULL, description = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits parameters of a NIG using maximum product spacings  

    # Example:
    #   set.seed(4711); x = rnig(500); mps = .nigFit.mps(x)@fit$estimate; mps
    
    # FUNCTION:

    # Transform:
    x.orig = x
    x = as.vector(x)
    if (scale) {
        SD = sd(x)
        x = x / SD }

    # Settings:
    CALL = match.call()

    # Parameter Estimation:
    obj <- function(x, y = x, trace) {
        if (abs(x[2]) >= x[1]) return(1e9)
        DH = diff(c(0, na.omit(.pnigC(sort(y), x[1], x[2], x[3], x[4])), 1))
        f = -mean(log(DH[DH > 0]))*length(y)
        # Print Iteration Path:
        if (trace) {
            cat("\n Objective Function Value:  ", -f)
            cat("\n Parameter Estimates:       ", x[1], x[2], x[3], x[4], "\n")
        }
        f }
    eps = 1e-10
    BIG = 1000
    r = nlminb(start = c(alpha, beta, delta, mu), objective = obj,
        lower = c(eps, -BIG, eps, -BIG), upper = BIG, y = x, trace = trace)
    names(r$par) <- c("alpha", "beta", "delta", "mu")
    
    # Standard Errors:
    hessian = tsHessian(x = r$par, fun = obj, y = x, trace = FALSE)
    colnames(hessian) = rownames(hessian) = names(r$par)
    varcov = solve(hessian)
    par.ses = sqrt(diag(varcov))
    if (scale) par.ses = par.ses / c(SD, SD, 1/SD, 1/SD)

    # Add Title and Description:
    if (is.null(title)) title = "NIG mps Parameter Estimation"
    if (is.null(description)) description = description()

    # Result:
    if (scale) {
        r$par = r$par / c(SD, SD, 1/SD, 1/SD)
        r$objective = obj(r$par, y = as.vector(x.orig), trace = trace)
    }
    fit = list(
        estimate = r$par, 
        minimum = -r$objective, 
        error = par.ses,
        code = r$convergence)

    # Optional Plot:
    if (doplot) {
        x = as.vector(x.orig)
        if (span == "auto") span = seq(min(x), max(x), length = 501)
        z = density(x, n = 100, ...)
        x = z$x[z$y > 0]
        y = z$y[z$y > 0]
        y.points = dnig(span, r$par[1], r$par[2], r$par[3], r$par[4])
        ylim = log(c(min(y.points), max(y.points)))
        if (add) {
            lines(x = span, y = log(y.points), col = "steelblue")
        } else {
            plot(x, log(y), xlim = c(span[1], span[length(span)]),
                ylim = ylim, type = "p", xlab = "x", ylab = "log f(x)", ...)
            title("NIG MPS Parameter Estimation")
            lines(x = span, y = log(y.points), col = "steelblue")
        }
    }

    # Return Value:
    new("fDISTFIT",
        call = as.call(CALL),
        model = "Normal Inverse Gaussian Distribution",
        data = as.data.frame(x.orig),
        fit = fit,
        title = as.character(title),
        description = description() )
}


# ------------------------------------------------------------------------------


.nigFit.vmps  <- 
function (x, alpha = 1, beta = 0, delta = 1, mu = 0,
    scale = TRUE, doplot = TRUE, add = FALSE, span = "auto", trace = TRUE, 
    title = NULL,description = NULL, ...)
{
    # A function implemented by Yohan Chalabi

    # Description:
    #   Fits parameters of a NIG using maximum product spacings  

    # Example:
    #   set.seed(4711); x = rnig(500); vmps = .nigFit.vmps(x)@fit$estimate; vmps
    
    # FUNCTION:

    # Transform:
    x.orig = x
    x = as.vector(x)
    if (scale) {
        SD = sd(x)
        x = x/SD
    }
    
    # Settings:
    CALL = match.call()
    
    # Parameter Estimation:
    obj <- function(x, y = x, trace) {
        if (abs(x[2]) >= x[1]) return(1e+9)
        DH = diff(c(0, na.omit(.pnigC(sort(y), x[1], x[2], x[3], x[4])), 1))
        f = log(var(DH[DH > 0]))    
        if (trace) {
            cat("\n Objective Function Value:  ", -f)
            cat("\n Parameter Estimates:       ", x[1], x[2], x[3], x[4], "\n")
        }
        f
    }
    eps = 1e-10
    BIG = 1000
    r = nlminb(
        start = c(alpha, beta, delta, mu), 
        objective = obj,
        lower = c(eps, -BIG, eps, -BIG), 
        upper = BIG, y = x,
        trace = trace)
    names(r$par) <- c("alpha", "beta", "delta", "mu")
    
    # Standard Errors:
    hessian = tsHessian(x = r$par, fun = obj, y = x, trace = FALSE)
    colnames(hessian) = rownames(hessian) = names(r$par)
    varcov = solve(hessian)
    par.ses = sqrt(diag(varcov))
    if (scale) par.ses = par.ses / c(SD, SD, 1/SD, 1/SD)
    
    # Add Title and Description:
    if (is.null(title)) title = "NIG varMPS Parameter Estimation"
    if (is.null(description)) description = description()
    
    # Result:
    if (scale) {
        r$par = r$par/c(SD, SD, 1/SD, 1/SD)
        r$objective = obj(r$par, y = as.vector(x.orig), trace = trace)
    }
    fit = list(
        estimate = r$par, 
        minimum = -r$objective, 
        error = par.ses,
        code = r$convergence)
    
    # Optional Plot:
    if (doplot) {
        x = as.vector(x.orig)
        if (span == "auto") span = seq(min(x), max(x), length = 501)
        z = density(x, n = 100, ...)
        x = z$x[z$y > 0]
        y = z$y[z$y > 0]
        y.points = dnig(span, r$par[1], r$par[2], r$par[3], r$par[4])
        ylim = log(c(min(y.points), max(y.points)))
        if (add) {
            lines(x = span, y = log(y.points), col = "steelblue")
        } else {
            plot(x, log(y), xlim = c(span[1], span[length(span)]),
                ylim = ylim, type = "p", xlab = "x", ylab = "log f(x)", ...)
            title("NIG varMPS Parameter Estimation")
            lines(x = span, y = log(y.points), col = "steelblue")
        }
    }
    
    
    # Return Value:
    new("fDISTFIT",
        call = as.call(CALL),
        model = "Normal Inverse Gaussian Distribution",
        data = as.data.frame(x.orig),
        fit = fit,
        title = as.character(title),
        description = description() )
}


################################################################################

