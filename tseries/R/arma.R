## Copyright (C) 1997-2000  Adrian Trapletti
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
## ARMA class
##

arma <-
function(x, order = c(1, 1), lag = NULL, coef = NULL,
         include.intercept = TRUE, series = NULL, qr.tol = 1e-07, ...)
{
    seqN <- function(N) {
        if(0==length(N)) NULL else if(N<=0) NULL else seq(N)
    }
  
    err <- function(coef) {
        u <- double(n)
        u[seqN(max.order)] <- 0
        u <- .C("arma",
                as.vector(x, mode = "double"),
                u = as.vector(u),
                as.vector(coef, mode = "double"),
                as.integer(lag$ar),
                as.integer(lag$ma),
                as.integer(ar.l),
                as.integer(ma.l),
                as.integer(max.order),
                as.integer(n),
                as.integer(include.intercept),
                PACKAGE="tseries")$u
        return(sum(u^2))
    }
  
    resid <- function(coef) {
        u <- double(n)
        u[seqN(max.order)] <- 0
        u <- .C("arma",
                as.vector(x, mode = "double"),
                u = as.vector(u),
                as.vector(coef, mode = "double"),
                as.integer(lag$ar),
                as.integer(lag$ma),
                as.integer(ar.l),
                as.integer(ma.l),
                as.integer(max.order),
                as.integer(n),
                as.integer(include.intercept),
                PACKAGE="tseries")$u
        return(u)
    }
  
    arma.init <- function() {
        k <- round(1.1*log(n))
        e <- na.omit(drop(ar.ols(x, order.max = k, aic = FALSE,
                                 demean = FALSE,
                                 intercept = include.intercept)$resid))
        ee <- embed(e, max.order+1)
        xx <- embed(x[-(1:k)], max.order+1)
        if(include.intercept == TRUE) {
            if(is.null(lag$ar))
                coef <- lm(xx[,1]~ee[,lag$ma+1])$coef
            else if(is.null(lag$ma))
                coef <- lm(xx[,1]~xx[,lag$ar+1])$coef
            else
                coef <- lm(xx[,1]~xx[,lag$ar+1]+ee[,lag$ma+1])$coef
            coef <- c(coef[-1], coef[1])
        } 
        else {
            if(is.null(lag$ar))
                coef <- lm(xx[,1]~ee[,lag$ma+1]-1)$coef
            else if(is.null(lag$ma))
                coef <- lm(xx[,1]~xx[,lag$ar+1]-1)$coef
            else
                coef <- lm(xx[,1]~xx[,lag$ar+1]+ee[,lag$ma+1]-1)$coef
        }
        return(coef) 
    }
  
    if(!is.null(order) && !is.null(lag))
        warning("order is ignored")
    if(is.null(order) && is.null(lag))
        stop("order or lag must be given")
    if(is.null(lag) && !is.null(order))
        lag <- list(ar=seqN(order[1]), ma=seqN(order[2]))
    lag$ar <- unique(lag$ar)
    lag$ma <- unique(lag$ma)
    max.order <- max(unlist(lag),0)
    ar.l <- length(lag$ar)
    ma.l <- length(lag$ma)
    if(NCOL(x) > 1)
        stop("x is not a vector or univariate time series")
    if(is.null(series)) series <- deparse(substitute(x))
    ists <- is.ts(x)
    x <- as.ts(x)
    xfreq <- frequency(x)
    if(any(is.na(x))) stop("NAs in x")
    if(ists) xtsp <- tsp(x)
    n <- length(x)
    if(!is.null(unlist(lag)))
        if((min(unlist(lag)) < 1) || (max(unlist(lag)) > (n-1)))
            stop("invalid lag")
    ncoef <- length(unlist(lag))+as.numeric(include.intercept)
    if(is.null(coef)) {
        if(!is.null(unlist(lag)))
            coef <- arma.init()
        else
            coef <- 0
    }
    if(length(coef) != ncoef) stop("invalid coef")
    md <- optim(coef, err, gr=NULL, hessian=TRUE, ...)
    coef <- md$par
    rank <- qr(md$hessian, qr.tol)$rank
    if(rank != ncoef) {
        vc <- matrix(NA, nrow = ncoef, ncol = ncoef)
        warning("singular Hessian")
    }
    else {
        vc <- 2*md$value/n*solve(md$hessian)
        if(any(diag(vc) < 0)) warning("Hessian negative-semidefinite")
    }
    e <- resid(coef)
    e[seqN(max.order)] <- NA
    f <- x-e
    if(ists) {
        attr(e, "tsp") <- xtsp
        attr(e, "class") <- "ts"
        attr(f, "tsp") <- xtsp
        attr(f, "class") <- "ts"
    }
    nam.ar <- if(!is.null(lag$ar))
        paste("ar", lag$ar, sep = "")
    else
        NULL
    nam.ma <- if(!is.null(lag$ma))
        paste("ma", lag$ma, sep = "")
    else
        NULL
    nam.int <- if(include.intercept) "intercept" else NULL
    nam.coef <- c(nam.ar, nam.ma, nam.int)
    names(coef) <- nam.coef
    colnames(vc) <- rownames(vc) <- nam.coef
    arma <- list(coef = coef,
                 css = md$value,
                 n.used = n,
                 residuals = e,
                 fitted.values = f,
                 series = series,
                 frequency = xfreq,
                 call = match.call(),
                 vcov = vc,
                 lag = lag,
                 convergence = md$convergence,
                 include.intercept = include.intercept)
    class(arma) <- "arma"
    return(arma)
}

coef.arma <-
function(object, ...)
{
    if(!inherits(object, "arma"))
        stop("method is only for arma objects")
    return(object$coef)
}

vcov.arma <-
function(object, ...)
{
    if(!inherits(object, "arma"))
        stop("method is only for arma objects")
    return(object$vcov)
}

residuals.arma <-
function(object, ...)
{
    if(!inherits(object, "arma"))
        stop("method is only for arma objects")
    return(object$residuals)
}

fitted.arma <-
function(object, ...)
{
    if(!inherits(object, "arma"))
        stop("method is only for arma objects")
    return(object$fitted.values)
}

print.arma <-
function(x, digits = max(3, getOption("digits") - 3), ...)
{
    if(!inherits(x, "arma"))
        stop("method is only for arma objects")
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("Coefficient(s):\n")
    print.default(format(coef(x), digits = digits),
                  print.gap = 2, quote = FALSE)
    cat("\n")
    invisible(x)
}

summary.arma <-
function(object, ...)
{
    if(!inherits(object, "arma"))
        stop("method is only for arma objects")
    ans <- NULL
    ans$residuals <- na.remove(object$residuals)
    tval <- object$coef / sqrt(diag(object$vcov))
    ans$coef <- cbind(object$coef, sqrt(diag(object$vcov)), tval,
                      2 * (1-pnorm(abs(tval))))
    dimnames(ans$coef) <-
        list(names(object$coef),
             c(" Estimate"," Std. Error"," t value","Pr(>|t|)"))
    ans$call <- object$call
    ans$nn <- object$nn
    ans$css <- object$css
    ans$var <- var(ans$residuals)
    ans$aic <- (object$n.used * (1+log(2*pi)) + object$n.used *
                log(ans$var) + 2 * length(object$coef))
    ans$p <- max(object$lag$ar, 0)
    ans$q <- max(object$lag$ma, 0)
    class(ans) <- "summary.arma"
    return(ans)
}

print.summary.arma <-
function(x, digits = max(3, getOption("digits") - 3),
         signif.stars = getOption("show.signif.stars"), ...)
{
    if(!inherits(x, "summary.arma"))
        stop("method is only for summary.arma objects")
    cat("\nCall:\n", deparse(x$call), "\n", sep = "")
    cat("\nModel:\nARMA(",x$p,",",x$q,")\n", sep = "")
    cat("\nResiduals:\n")
    rq <- structure(quantile(x$residuals),
                    names = c("Min","1Q","Median","3Q","Max"))
    print(rq, digits = digits, ...)
    cat("\nCoefficient(s):\n")
    printCoefmat(x$coef, digits = digits,
                 signif.stars = signif.stars, ...)
    cat("\nFit:\n")
    cat("sigma^2 estimated as ", format(x$var, digits = digits), 
        ",  Conditional Sum-of-Squares = ", format(round(x$css, 2)), 
        ",  AIC = ", format(round(x$aic, 2)), "\n", sep = "")
    cat("\n")
    invisible(x)
}

plot.arma <-
function(x, ask = interactive(), ...)
{
    if(!inherits(x, "arma"))
        stop("method is only for arma objects")
    op <- par()
    par(ask = ask, mfrow = c(2, 1))
    data <- eval.parent(parse(text = x$series))
    if(any(is.na(data))) stop(paste("NAs in", x$series))
    plot(data, main = x$series, ylab = "Series")
    plot(x$residuals, main = "Residuals", ylab = "Series")
    acf(data, main = paste("ACF of", x$series))
    acf(x$residuals, main = "ACF of Residuals",
        na.action = na.remove)
    pacf(data, main = paste("PACF of", x$series))
    pacf(x$residuals, main = "PACF of Residuals",
         na.action = na.remove)
    par(ask = op$ask, mfrow = op$mfrow)
    invisible(x)
}
