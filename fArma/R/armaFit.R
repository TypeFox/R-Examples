
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
# FUNCTION:               PARAMETER ESTIMATION:
#  armaFit                 Fits parameters for ARMA Time Series process
#  .arFit                  Internal function called by armaFit
#  .arimaFit               Internal function called by armaFit
#  .arfimaFit              Internal function called by armaFit
################################################################################


armaFit =
function(
formula, data, method = c("mle", "ols"),
include.mean = TRUE, fixed = NULL, title = NULL, description = NULL, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits Model Parameters for an ARMA Time Series Process

    # Arguments:
    #   method - "mle", "ols"

    # Notes:
    #   Valid formulas are:
    #       "ar" "ma" "arma" "arima" "arfima" (not documented "fracdiff")
    #       Note, "arma(p,q)" uses arima(p,0,q)

    # Details:
    #   R-base:
    #       arima(
    #           x,
    #           order = c(0, 0, 0),
    #           seasonal = list(order = c(0, 0, 0), period = NA),
    #           xreg = NULL,
    #           include.mean = TRUE,
    #           transform.pars = TRUE,
    #           fixed = NULL,
    #           init = NULL,
    #           method = c("CSS-ML", "ML", "CSS"),
    #           n.cond,
    #           optim.control = list(),
    #           kappa = 1e+06)
    #   Compare with SPlus:
    #       arima.mle(
    #           x,
    #           model,
    #           n.cond,
    #           xreg=NULL,
    #           ...)

    # FUNCTION:

    # Check for Method:
    # ar.method       = c("yw", "burg1", "burg2", "ols", "mle")
    # arma.method     = c("CSS")
    # arima.method    = c("CSS-ML", "ML", "CSS")
    # fracdiff.method = NA
    method = method[1]  # Don't use match.arg(methods)

    # Call:
    fit = NULL
    call = match.call()

    # Add to Fracdiff: h and M default Settings
    mf = match.call(expand.dots = TRUE)
    m = match("h", names(mf), 0)
    if (m == 0) h = -1 else h = eval(mf[[m]])
    fracdiff.h = h
    m = match("M", names(mf), 0)
    if (m == 0) M = 100 else M = eval(mf[[m]])
    fracdiff.M = M

    # Get Series:
    # DW 2005-09-03
    mf = match.call(expand.dots = FALSE)
    m = match(c("formula", "data"), names(mf), 0)
    mf = mf[c(1, m)]
    mf[[1]] = as.name(".modelSeries")
    mf$fake = FALSE
    mf$lhs = TRUE
    if (missing(data)) data = eval(parse(text = search()[2]), parent.frame())
    mf$data = as.matrix(data)
    x = eval(mf, parent.frame())
    x = ts = as.vector(x[, 1])

###     # Allow for univariate 'timeSeries' Objects:
###     # Added 2004-09-04 DW
###     if (class(ts) == "timeSeries") ts = as.vector(ts)

    # Which Model?
    # DW 2006-02-21
    K = length(formula)
    end = regexpr("\\(", as.character(formula[K]))-1
    tsmodel =  substr(as.character(formula[K]), 1, end)

    # Valid Model?
    valid = FALSE
    if (tsmodel == "ar" ) valid = TRUE
    if (tsmodel == "ma" ) valid = TRUE
    if (tsmodel == "arma") valid = TRUE
    if (tsmodel == "arima") valid = TRUE
    if (tsmodel == "arfima") valid = TRUE
    if (tsmodel == "fracdiff") valid = TRUE
    if (!valid) stop("Invalid Formula Specification")

    # Which Order?
    start = regexpr("\\(", as.character(formula[K]))+1
    end   = regexpr("\\)", as.character(formula[K]))-1
    order = substr(as.character(formula[K]), start, end)

    if (tsmodel == "arfima" || tsmodel == "fracdiff") {
        # method will be ignored ...
        pos = regexpr(",", order)
        p = as.integer(substr(order, 1, pos-1))
        q = as.integer(substr(order, pos+1, nchar(order)))
        order = c(p, q)
        tsmodel = "arfima"
    }

    if (tsmodel == "arima") {
        if (method == "mle") method = "CSS-ML"
        pos = regexpr(",", order)
        p = as.integer(substr(order, 1, pos-1))
        order = substr(order, pos+2, nchar(order))
        d = as.integer(substr(order, 1, pos-1))
        q = as.integer(substr(order, pos+1, nchar(order)))
        order = c(p = p, d = d, q = q)
    }

    if (tsmodel == "arma") {
        if (method == "mle") method = "CSS-ML"
        # "arma" uses "arima"
        pos = regexpr(",", order)
        p = as.integer(substr(order, 1, pos-1))
        q = as.integer(substr(order, pos+1, nchar(order)))
        order = c(p = p, d = 0, q = q)
        tsmodel = "arima"
    }

    if (tsmodel == "ar") {
        # if method is CSS-ML, CSS, or ML, then "ar" uses "arima":
        order = as.integer(order)
        if (method == "mle") method = "CSS-ML"
        if (method == "CSS-ML" | method == "CSS" | method == "ML") {
            p = order
            order = c(p = p , d = 0, q = 0)
            tsmodel = "arima"
        }
    }

    if (tsmodel == "ma") {
        # if method is CSS-ML, CSS, or ML, then "ma" uses "arima":
        if (method == "mle") method = "CSS-ML"
        order = as.integer(order)
        order = c(p = 0 , d = 0, q = order)
        tsmodel = "arima"
    }

    # Which Function?
    fun = match.fun(paste(".", tsmodel, "Fit", sep = ""))

    # Fit:
    fit = fun(x = ts, order = order, include.mean = include.mean,
        method = method[1], fixed = fixed, ...)

    # "ols" specific:
    if (method == "ols") {
        se.coef = unlist(fit$se.coef)
        if (include.mean){
            ols.mean = se.coef[1]
            fit$se.coef = c(se.coef[-1], ols.mean) }
    }
    fit$call = call
    fit$tsmodel = tsmodel
    fit$class = "fARMA"
    class(fit) = "list"

    # Add title and desription:
    if (is.null(title)) title = "ARIMA Modelling"
    if (is.null(description)) description = description()


    # Parameters:
    parameter = list(include.mean = include.mean, fixed = fixed)
    if (tsmodel == "arfima" || tsmodel == "fracdiff") {
        parameter$M = fracdiff.M
        parameter$h = fracdiff.h
    }

    # Return Value:
    new("fARMA",
        call = as.call(match.call()),
        formula = as.formula(formula),
        method = as.character(method),
        parameter = parameter,
        data = list(x = x),
        fit = fit,
        residuals = list(residuals = fit$residuals),
        fitted = list(fitted = fit$fitted.values),
        title = as.character(title),
        description = as.character(description) )
}


# ------------------------------------------------------------------------------


.arFit =
function(x, order, include.mean, fixed = NULL,
method = c("yw", "burg1", "burg2", "ols", "mle"), M = NULL, h = NULL, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits an AR time series model

    # Note:
    #   Calls ar() from R-stats.

    # FUNCTION:

    # Fit:
    call = match.call()
    var.method = as.integer(substring(method[1], 5, 5))
    method = substr(method[1], 1, 4)
    fit = ar(x = x, aic = FALSE, order.max = order, method = method,
        var.method = var.method, demean = include.mean,
        intercept = include.mean, ...)

    # Add and Modify:
    fit$call = call
    fit$tstitle = paste("AR(",
        as.character(order), ") with method: ", method, sep = "")
    fit$order = order

    # Residuals:
    fit$residuals = fit$resid
    fit$fitted.values = x - fit$resid
    fit$sigma2 = fit$var.pred

    # Coefficients:
    if (method == "ols") {
        fit$coef = fit$ar[,,1]
    } else {
        fit$coef = fit$ar
    }
    names(fit$coef) = c(paste("ar", 1:order, sep=""))

    # Mean:
    if (include.mean) {
        coeff = c(fit$coef, fit$x.mean)
        names(coeff) = c(names(fit$coef), "intercept")
        fit$coef = coeff
    }

    # Standard Errors:
    if (method == "ols") {
        fit$se.coef = fit$asy.se.coef
        n = sqrt(length(as.vector(fit$se.coef)))
        fit$var.coef = matrix(rep(NA, times = n*n), ncol = n)
    } else {
        fit$var.coef = fit$asy.var.coef
        fit$se.coef = sqrt(diag(fit$asy.var.coef))
        if (include.mean) {
            m = dim(fit$asy.var.coef)[1] + 1
            var.coef = matrix(rep(NA, times = m*m), m, m)
            for ( i in 1:(m-1) ) {
                for ( j in 1:(m-1) ) {
                    var.coef[i,j] = fit$var.coef[i,j]
                }
            }
            fit$var.coef = var.coef
            fit$se.coef = c(fit$se.coef, NA)
        }
    }

    # Add Data:
    fit$x = x

    # Return Value:
    fit
}


# ------------------------------------------------------------------------------


.arimaFit =
function (x, order, include.mean, fixed,
method = c("CSS-ML", "ML", "CSS"), M = NULL, h = NULL, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits an ARIMA time series model

    # Note:
    #   Calls arima() from R-stats.

    # FUNCTION:

    # Fit:
    call = match.call()
    fit = arima(x = x, order = order, method =  method[1],
        include.mean = include.mean, fixed = fixed, ...)

    # Add Title:
    fit$tstitle = paste("ARIMA(",
        as.character(order[1]), ",", as.character(order[2]), ",",
        as.character(order[3]), ") with method: ", method[1], sep = "")

    # Add Data:
    fit$x = x

    # Add Fitted Values:
    fit$fitted.values = fit$x - fit$residuals

    # Add Standard Errors:
    fit$se.coef = sqrt(diag(fit$var.coef))

    # Add Call:
    fit$call = call

    # Return Value:
    fit
}


# ------------------------------------------------------------------------------


.arfimaFit =
function (x, order, include.mean, fixed, method = "arfima",
M = 100, h = -1, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits an ARFIMA (FRACDIFF) time series model

    # Arguments:
    #   x - time series for the ARIMA model
    #   nar - number of autoregressive parameters
    #   nma - number of moving average parameters
    #   ar - initial autoregressive parameters
    #   ma - initial moving average parameters
    #   dtol - desired accurcay for d, by default (and if
    #       negative), (4th root of machine precision)
    #       is used.  dtol will be changed internally if
    #       necessary
    #   drange - interval over which the likelihood function is
    #       to be maximized as a function of d
    #   h - finite difference interval
    #   M - number of terms in the likelihood approximation
    #       (see Haslett and Raftery 1989)

    # Note:
    #   A Builtin Copy from R's fracdiff Package
    #   Calls fracdiff() from R-fracdiff

    # FUNCTION:

    # Settings:
    call = match.call()
    nar = order[1]
    nma = order[2]
    ar = rep(NA, max(order[1], 1))
    ma = rep(NA, max(order[2], 1))
    dtol = .Machine$double.eps^0.25 # ~ 1.22e-4
    drange = c(0, 0.5)

    # fracdiff:
    if (any(is.na(x)))
        stop("missing values not allowed in time series")
    if (is.matrix(x) && ncol(x) > 2)
        stop("multivariate time series not allowed")
    n = length(x)
    npq = nar + nma
    npq1 = npq + 1
    lwork = max(npq+2*(n+M), 3*n+(n+6)*npq+npq%/%2+1, (3+2*npq1)*npq1+1)
    ar[is.na(ar)] = 0
    ma[is.na(ma)] = 0

    # if dtol < 0: the fortran code will choose defaults
    result = .Fortran("fracdf", as.double(x), as.integer(n),
        as.integer(M), as.integer(nar), as.integer(nma),
        dtol = as.double(dtol), drange = as.double(drange),
        hood = double(1), d = double(1), ar = as.double(ar),
        ma = as.double(ma), w = double(lwork), as.integer(lwork),
        info = integer(1), .Machine$double.xmin,
        .Machine$double.xmax, .Machine$double.neg.eps,
        .Machine$double.eps, PACKAGE = "fArma")
    if (result$info) switch(result$info,
        stop("insufficient workspace"),
        stop("error in gamma function"),
        stop("invalid MINPACK input"),
        warning(" Warning in gamma function"),
        warning(" Optimization failure"),
        warning(" Optimization limit reached"))
    hess = .Fortran("fdhpq",
         as.double(x), hess = double(npq1 * npq1), as.integer(npq1),
         result$w, PACKAGE = "fArma")$hess
    temp = .Fortran("fdcov", as.double(x), as.double(result$d),
         h = as.double(if (missing(h)) -1 else h), hd = double(npq1),
         cov = hess, as.integer(npq1), cor = hess, as.integer(npq1),
         se = double(npq1), result$w, info = integer(1),
         PACKAGE = "fArma")
    if (temp$info) switch(temp$info,
         warning(" Warning in gamma function"),
         warning(" Singular Hessian matrix"),
         warning(" Unable to compute correlation matrix"),
         stop("error in gamma function"))
    if (npq == 0) {
        result$ar = NULL
        result$ma = NULL }
    nam = "d"
    if (nar) nam = c(nam, paste("ar", 1:nar, sep = ""))
    if (nma) nam = c(nam, paste("ma", 1:nma, sep = ""))
    hess = matrix(hess, nrow = npq1, ncol = npq1, dimnames = list(nam, nam))
    hess[1, ] = temp$hd
    hess[row(hess) > col(hess)] = hess[row(hess) < col(hess)]
    se.ok = temp$info != 0 || temp$info < 3

    # Fitting Result:
    fit = list(
        log.likelihood = result$hood,
        d = result$d,
        ar = result$ar, ma = result$ma,
        covariance.dpq = array(temp$cov, c(npq1, npq1), list(nam, nam)),
        stderror.dpq = if (se.ok) temp$se, # else NULL
        correlation.dpq = if (se.ok) array(temp$cor, c(npq1, npq1)), # else NULL
        h = temp$h, d.tol = result$dtol, M = M, hessian.dpq = hess)

    # Add ts Title:
    fit$tstitle = paste("FRACDIFF(", as.character(order[1]), ",",
        as.character(order[2]), ") with method: ", method[1], sep = "")

    # Add Series:
    fit$x = x

    # Add Coefficients:
    fit$coef = c(fit$d, fit$ar, fit$ma)
    namesCoef = "d"
    if (order[1] > 0) {
        names.ar = c(paste("ar", 1:order[1], sep=""))
        namesCoef = c(namesCoef, names.ar) }
    if (order[2] > 0) {
        names.ma = c(paste("ma", 1:order[2], sep=""))
        namesCoef = c(namesCoef, names.ma) }
    names(fit$coef) = namesCoef
    fit$var.coef = fit$correlation.dpq

    # Add Fitted Values:
    n = 0:fit$M
    w = lgamma(-fit$d+n) - (lgamma(-fit$d)+lgamma(n+1))
    w = exp(w)
    fit$fitted.values = filter(fit$x, w, sides = 1)

    # Add Residuals:
    fit$residuals = x - fit$fitted.values

    # Add Standard Errors:
    fit$se.coef = fit$stderror.dpq

    # Add fracdiff Parameters:
    fit$fracdiff = c(M, h)

    # Add Call:
    fit$call = call

    # Return Value:
    fit
}


################################################################################

