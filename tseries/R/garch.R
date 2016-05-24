## Copyright (C) 1997-1999  Adrian Trapletti
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## A copy of the GNU General Public License is available via WWW at
## http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
## writing to the Free Software Foundation, Inc., 59 Temple Place,
## Suite 330, Boston, MA  02111-1307  USA.

##
## GARCH class
##

garch <-
function (x, order = c(1, 1), series = NULL, control = garch.control(...), ...)
{
    if(NCOL(x) > 1)
        stop("x is not a vector or univariate time series")
    if(!is.vector(order)) stop("order is not a vector")
    switch(control$grad,
           analytical = (agrad <- TRUE),
           numerical = (agrad <- FALSE))
    if(is.null(series)) series <- deparse(substitute(x))
    ists <- is.ts(x)
    x <- as.ts(x)
    xfreq <- frequency(x)
    if(any(is.na(x))) stop("NAs in x")
    if(ists) xtsp <- tsp(x)
    x <- as.matrix(x)
    n <- nrow(x)
    e <- double(n)
    ncoef <- order[1]+order[2]+1
    hess <- matrix(0.0, ncoef, ncoef)
    small <- 0.05
    coef <- control$start
    if(is.null(coef))
        coef <- c(var(x)*(1.0-small*(ncoef-1)),rep.int(small,ncoef-1))
    if(!is.vector(coef)) stop("coef is not a vector")
    if(ncoef != length(coef)) stop("incorrect length of coef")
    nlikeli <- 1.0e+10
    fit <- .C("fit_garch",
              as.vector(x, mode = "double"),
              as.integer(n),
              coef = as.vector(coef, mode = "double"),
              as.integer(order[1]),
              as.integer(order[2]),
              as.integer(control$maxiter),
	      as.double(control$abstol),
	      as.double(control$reltol),
	      as.double(control$xtol),
	      as.double(control$falsetol),
              nlikeli = as.double(nlikeli),
              as.integer(agrad),
              as.integer(control$trace),
              PACKAGE="tseries")
    pred <- .C("pred_garch",
               as.vector(x, mode = "double"),
               e = as.vector(e, mode = "double"),
               as.integer(n),
               as.vector(fit$coef, mode = "double"),
               as.integer(order[1]),
               as.integer(order[2]),
               as.integer(FALSE),
               PACKAGE = "tseries")
    com.hess <- .C("ophess_garch",
                   as.vector(x, mode = "double"),
                   as.integer(n),
                   as.vector(fit$coef, mode = "double"),
                   hess = as.matrix(hess),
                   as.integer(order[1]),
                   as.integer(order[2]),
                   PACKAGE="tseries")
    rank <- do.call("qr", c(list(x = com.hess$hess), control$qr))$rank
    if(rank != ncoef) {
	vc <- matrix(NA, nrow = ncoef, ncol = ncoef)
        warning("singular information")
    }
    else
        vc <- solve(com.hess$hess)
    sigt <- sqrt(pred$e)
    sigt[1:max(order[1],order[2])] <- rep.int(NA, max(order[1],order[2]))
    f <- cbind(sigt,-sigt)
    colnames(f) <- c("sigt","-sigt")
    e <- as.vector(x)/sigt  
    if(ists) {
        attr(e, "tsp") <-  attr(f, "tsp") <- xtsp
        attr(e, "class") <- attr(f, "class") <- "ts"
    }
    names(order) <- c("p","q")
    coef <- fit$coef
    nam.coef <- "a0"
    if(order[2] > 0)
        nam.coef <- c(nam.coef, paste("a", seq(order[2]), sep = ""))
    if(order[1] > 0)
        nam.coef <- c(nam.coef, paste("b", seq(order[1]), sep = ""))
    names(coef) <- nam.coef
    colnames(vc) <- rownames(vc) <- nam.coef
    garch <- list(order = order,
                  coef = coef,
                  n.likeli = fit$nlikeli,
                  n.used = n,
                  residuals = e,
                  fitted.values = f,
                  series = series, 
                  frequency = xfreq,
                  call = match.call(),
		  vcov = vc)
    class(garch) <- "garch"
    return(garch)
}

garch.control <-
function(maxiter = 200, trace = TRUE, start = NULL, grad = c("analytical","numerical"),
  abstol = max(1e-20, .Machine$double.eps^2),
  reltol = max(1e-10, .Machine$double.eps^(2/3)),
  xtol = sqrt(.Machine$double.eps),
  falsetol = 1e2 * .Machine$double.eps, ...)
{
  rval <- list(maxiter = maxiter, trace = trace, start = start, grad = match.arg(grad),
    abstol = abstol, reltol = reltol, xtol = xtol, falsetol = falsetol)
  rval$qr <- list(...)
  rval
}

coef.garch <-
function(object, ...)
{
    if(!inherits(object, "garch"))
        stop("method is only for garch objects")
    return(object$coef)
}

vcov.garch <-
function(object, ...)
{
    if(!inherits(object, "garch"))
        stop("method is only for garch objects")
    return(object$vcov)
}

residuals.garch <-
function(object, ...)
{
    if(!inherits(object, "garch"))
        stop("method is only for garch objects")
    return(object$residuals)
}

fitted.garch <-
function(object, ...)
{
    if(!inherits(object, "garch"))
        stop("method is only for garch objects")
    return(object$fitted.values)
}

print.garch <-
function(x, digits = max(3, getOption("digits") - 3), ...)
{
    if(!inherits(x, "garch"))
        stop("method is only for garch objects")
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("Coefficient(s):\n")
    print.default(format(coef(x), digits = digits), print.gap = 2,
                  quote = FALSE)
    cat("\n")
    invisible(x)
}

summary.garch <-
function(object, ...)
{
    if(!inherits(object, "garch"))
        stop("method is only for garch objects")
    ans <- NULL
    ans$residuals <- na.remove(object$residuals)
    tval <- object$coef / sqrt(diag(object$vcov))
    ans$coef <- cbind(object$coef, sqrt(diag(object$vcov)), tval,
                      2*(1-pnorm(abs(tval))))
    dimnames(ans$coef) <-
        list(names(object$coef),
             c(" Estimate"," Std. Error"," t value","Pr(>|t|)"))
    ans$call <- object$call
    ans$order <- object$order
    Residuals <- ans$residuals
    ans$j.b.test <- jarque.bera.test(Residuals)
    Squared.Residuals <- ans$residuals^2
    ans$l.b.test <- Box.test(Squared.Residuals, type = "Ljung-Box")
    class(ans) <- "summary.garch"
    return(ans)
}

plot.garch <- function(x, ask = interactive(), ...)
{
    if(!inherits(x, "garch"))
        stop("method is only for garch objects")
    op <- par()
    par(ask = ask, mfrow = c(2,1))
    data <- eval.parent(parse(text=x$series))
    if(any(is.na(data))) stop(paste("NAs in", x$series))
    plot(data, main = x$series, ylab = "Series")
    plot(x$residuals, main = "Residuals", ylab = "Series")
    hist(data, main = paste("Histogram of", x$series), xlab = "Series")
    hist(x$residuals, main = "Histogram of Residuals", xlab = "Series")
    qqnorm(data, main = paste("Q-Q Plot of", x$series),
           xlab = "Normal Quantiles")
    qqnorm(x$residuals, main = "Q-Q Plot of Residuals",
           xlab = "Normal Quantiles")
    acf(data^2, main = paste("ACF of Squared", x$series))
    acf(x$residuals^2, main = "ACF of Squared Residuals",
        na.action = na.remove)
    par(ask = op$ask, mfrow = op$mfrow)
    invisible(x)
}

print.summary.garch <-
function(x, digits = max(3, getOption("digits") - 3),
         signif.stars = getOption("show.signif.stars"), ...)
{
    if(!inherits(x, "summary.garch"))
        stop("method is only for summary.garch objects")
    cat("\nCall:\n", deparse(x$call), "\n", sep = "")
    cat("\nModel:\nGARCH(", x$order[1], ",", x$order[2], ")", "\n",
        sep = "")
    cat("\nResiduals:\n")
    rq <- structure(quantile(x$residuals),
                    names = c("Min","1Q","Median","3Q","Max"))
    print(rq, digits = digits, ...)
    cat("\nCoefficient(s):\n")
    printCoefmat(x$coef, digits = digits,
                 signif.stars = signif.stars, ...)
    cat("\nDiagnostic Tests:")
    print(x$j.b.test)
    print(x$l.b.test)
    invisible(x)
}

predict.garch <-
function(object, newdata, genuine = FALSE, ...)
{
    if(!inherits(object, "garch"))
        stop("method is only for garch objects")
    if(missing(newdata)) {
        newdata <- eval.parent(parse(text=object$series))
        if(any(is.na(newdata))) stop("NAs in newdata")
    }
    if(NCOL(newdata) > 1)
        stop("newdata is not a vector or univariate time series")
    ists <- is.ts(newdata)
    if(ists) newdata.tsp <- tsp(newdata)
    newdata <- as.matrix(newdata)
    n <- nrow(newdata)
    if(genuine) h <- double(n+1)
    else h <- double(n)
    pred <- .C("pred_garch",
               as.vector(newdata, mode = "double"),
               h = as.vector(h, mode = "double"),
               as.integer(n),
               as.vector(object$coef, mode = "double"),
               as.integer(object$order[1]),
               as.integer(object$order[2]),
               as.integer(genuine),
               PACKAGE="tseries")
    pred$h <- sqrt(pred$h)
    pred$h[1:max(object$order[1],object$order[2])] <-
        rep.int(NA, max(object$order[1],object$order[2]))
    pred$h <- cbind(pred$h,-pred$h)
    if(ists) {
        attr(pred$h, "tsp") <-
            if(genuine)
                c(newdata.tsp[1],
                  newdata.tsp[2] + 1 / newdata.tsp[3],
                  newdata.tsp[3])
            else
                newdata.tsp
        attr(pred$h, "class") <- "ts"
    }
    return(pred$h)
}

logLik.garch <-
function(object, ...)
{
    if(!inherits(object, "garch"))
        stop("method is only for garch objects")
    n <- length(na.remove(object$residuals))
    val <- (-object$n.likeli) - 0.5*n*log(2*pi)
    attr(val, "df") <- length(object$coef)
    class(val) <- "logLik"
    return(val)
}
