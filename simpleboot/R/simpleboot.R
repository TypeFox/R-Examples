## Functions for simple bootstrapping procedures
##
## Copyright 2004--2007 Roger D. Peng <rpeng@jhsph.edu>

## (6/18/2003) The `two.boot' function was rewritten to (properly) use
## the strata argument of `boot'.  This simplifies the two sample
## bootstrap and was what I wanted to do all along. -RDP

## (5/27/03) This package is poorly designed because, among other
## things, it relies heavily on knowledge of the underlying
## representation of the class `boot'.  Also, the print and summary
## methods are kind of hacked. -RDP


hist.simpleboot <- function(x, do.rug = FALSE, xlab = "Bootstrap samples",
                            main = "", ...) {
    hist(x$t[, 1], xlab = xlab, main = main, ...)
    if(do.rug)
        rug(x$t[, 1])
    abline(v = ifelse(is.matrix(x$t0), x$t0[,1], x$t0[1]), col = 2, lty=3)

    invisible()
}


########################################################################

## Utility functions

## Extract bootstrap samples from the boot object.  Is this useful?
## 'idx' indicates which sample to extract.  Can be a vector.

perc.lm <- function(lm.boot.obj, p) {
    cmat <- samples(lm.boot.obj, "coef")
    apply(cmat, 1, quantile, probs = p)
}

perc <- function(boot.out, p = c(0.025, 0.975)) {
    if(inherits(boot.out, "lm.simpleboot"))
        return( perc.lm(boot.out, p) )
    if(!inherits(boot.out, "simpleboot"))
        stop("only use this function on 'simpleboot' objects")
    if(any(p < 0) || any(p > 1))
        stop("probabilities in 'p' must be between 0 and 1")
    if(!is.null(boot.out$student) && boot.out$student)
        x <- boot.out$t[, 1, drop = FALSE]
    else
        x <- boot.out$t
    out <- drop(apply(x, 2, quantile, probs = p))
    out
}


########################################################################

## Actual bootstrapping functions
## If 'quantile' is bootstrapped, then the 'probs' argument must be set

one.boot <- function(data, FUN, R, student = FALSE, M, weights = NULL, ...) {
    func.name <- ifelse(is.character(FUN), FUN, deparse(substitute(FUN)))
    extra <- list(...)
    
    if(func.name == "quantile") {
        if(is.na(match("probs", names(extra))))
            stop("'probs' argument must be specified")
        if(length(extra$probs) > 1)
            stop("can only bootstrap a single quantile")
    }
    func <- match.fun(FUN)

    boot.func <- function(x, idx) {
        fval <- func(x[idx], ...)

        if(student) {
            rs.x <- x[idx]
            b <- one.boot(rs.x, FUN, R = M, student = FALSE, M = NULL, weights = NULL, ...)
            fval <- c(fval, var(b$t))
        }
        fval
    }
    b <- boot(data, statistic = boot.func, R = R, weights = weights)
    b$student <- student
    structure(b, class = "simpleboot")
}



## FUN should be some sort of function which takes one argument, like
## 'mean' or 'median'.  Same as for the one sample bootstrap.

two.boot <- function(sample1, sample2, FUN, R, student = FALSE, M,
                     weights = NULL, ...) {
    func.name <- ifelse(is.character(FUN), FUN, deparse(substitute(FUN)))
    func <- match.fun(FUN)
    ind <- c(rep(1, length(sample1)), rep(2, length(sample2)))
    nobsgrp <- as.numeric(table(ind))
    extra <- list(...)

    if(func.name == "quantile") {
        if(is.na(match("probs", names(extra))))
            stop("'probs' argument must be specified")
        if(length(extra$probs) > 1)
            stop("can only bootstrap a single quantile")
    }
    boot.func <- function(x, idx) {
        d1 <- x[idx[ind == 1]]
        d2 <- x[idx[ind == 2]]
        fval <- func(d1, ...) - func(d2, ...)

        if(student) {
            b <- two.boot(d1, d2, FUN, R = M, student = FALSE, M = NULL,
                          weights = NULL, ...) 
            fval <- c(fval, var(b$t))
        }
        fval
    }
    if(!is.null(weights))
        weights <- unlist(weights)
    b <- boot(c(sample1, sample2), statistic = boot.func, R = R,
              weights = weights, strata = ind)
    b$student <- student
    b$student <- student
    structure(b, class = "simpleboot")
}


## x and y should either be vectors, or x can be a two-column matrix
## There are two options for FUN.
## (a) FUN takes to arguments:  FUN <- function(x, y) { ... }
## (b) FUN takes one argument that's either a two-column matrix or a
##     two-column data frame:
##     FUN <- function(x) { ... }

pairs.boot <- function(x, y = NULL, FUN, R, student = FALSE, M,
                       weights = NULL, ...) {
    func <- match.fun(FUN)
    
    if(is.null(y)) {
        if(!is.matrix(x) && !is.data.frame(x))
            stop("'x' must be a matrix or a data frame")
        data <- x
    }
    else {
        if(length(x) != length(y))
            stop("length of 'x' must equal length of 'y'")
        data <- cbind(x, y)
    }
    if(student && missing(M))
        stop("need to specify 'M' for studentized bootstrap")
    
    boot.func <- function(x, idx) {
        rs.x <- x[idx, ]

        if(is.null(y)) 
            fval <- func(rs.x, ...)        
        else 
            fval <- func(rs.x[,1], rs.x[,2], ...)
        if(student) {
            b <- pairs.boot(rs.x[,1], rs.x[,2], FUN, M, student = FALSE,
                            M = NULL, weights = NULL, ...)
            fval <- c(fval, var(b$t))
        }
        fval
    }
    b <- boot(data, statistic = boot.func, R = R, weights = weights)
    b$student <- student
    structure(b, class = "simpleboot")
}














summary.lm.simpleboot <- function(object, ...) {
    summary.object <- object
    class(summary.object) <- "summary.lm.simpleboot"
    params <- sapply(object$boot.list, "[[", "coef")        
    summary.object$stdev.params <- apply(params, 1, sd)

    summary.object
}

print.summary.lm.simpleboot <- function(x, ...) {    
    print.lm.simpleboot(x)
    
    cat("Bootstrap SD's:\n")
    print.default(format(x$stdev.params), print.gap = 2, quote = FALSE)
    cat("\n")
}

print.lm.simpleboot <- function(x, ...) {
    cat("BOOTSTRAP OF LINEAR MODEL  (method = ", x$method, ")\n\n", sep = "")
    cat("Original Model Fit\n")
    cat("------------------")
    print(x$orig.lm)
}

fitted.lm.simpleboot <- function(object, ...) {
    samples(object, "fitted")
}

model.frame.lm.simpleboot <- function(formula, ...) {
    model.frame(formula$orig.lm)
}

plot.lm.simpleboot <- function(x, add = FALSE, ...) {
    if(ncol(model.frame(x)) > 2)
        stop("cannot plot bootstrap regressions with dimension > 2")
    xpts <- x$new.xpts    
    bmat <- sapply(x$boot.list, "[[", "fitted")
    std <- apply(bmat, 1, sd, na.rm = T)
    mask.response <- -attr(terms(x$orig.lm), "response")
    orig.pred <- predict(x$orig.lm, xpts)
    
    if(!add) {
        mf <- model.frame(x)        
        xdata <- mf[, mask.response]
        ydata <- mf[, -mask.response]               
        plot(xdata, ydata, ...)
    }
    abline(x$orig.lm)
    lines(xpts[,1], orig.pred + 2*std, lty=3)
    lines(xpts[,1], orig.pred - 2*std, lty=3)

    invisible()          
}

samples <- function(object, name = c("fitted", "coef", "rsquare", "rss")) {
    name <- match.arg(name)
    
    if(!inherits(object, c("lm.simpleboot", "loess.simpleboot")))
        stop("only use with 'lm.simpleboot' or 'loess.simpleboot' object")
    boot.list <- object$boot.list
    
    if(is.null(boot.list[[1]][[name]]))
        stop(gettextf("bootstrap model does not have '%s' values", name))
    sapply(boot.list, "[[", name)
}

lm.boot <- function(lm.object, R, rows = TRUE, new.xpts = NULL, ngrid = 100,
                    weights = NULL) {
    orig.data <- model.frame(lm.object)

    if(ncol(orig.data) == 2 && is.null(new.xpts)) {
        mask.response <- -attr(terms(lm.object), "response")
        range.x <- range(orig.data[, mask.response])
        new.xpts <- data.frame(seq(range.x[1], range.x[2], len = ngrid))
        names(new.xpts) <- names(orig.data)[mask.response]
    }
    if(is.null(weights))
        weights <- rep(1, NROW(orig.data))
    boot.list <- lm.boot.resample(lm.object, R, rows, new.xpts, weights)
    
    structure(list(method = ifelse(rows, "rows", "residuals"),
                   boot.list = boot.list, orig.lm = lm.object,
                   new.xpts = new.xpts, weights = weights),
              class = "lm.simpleboot")
}

lm.boot.resample <- function(lm.object, R, rows, new.xpts, weights) {
    boot.list <- vector("list", length = R)
    names(boot.list) <- as.character(1:R)
    yhat <- predict(lm.object)
    f <- formula(lm.object)
    mframe <- model.frame(lm.object)
    nobs <- nrow(mframe)

    for(i in 1:R) {
        if(rows) {
            boot.idx <- sample(1:nobs, replace = TRUE, prob = weights)
            mf <- mframe[boot.idx, ]
        }
        else {
            mf <- model.frame(lm.object)
            rstar <- sample(residuals(lm.object), replace = TRUE,
                            prob = weights)

            ## Generate new responses with resampled residuals
            mf[[attr(terms(lm.object), "response")]] <- yhat + rstar
        }
        ## rs.lm <- lm(f, data = mf)
        rs.lm <- update(lm.object, data = mf)
        
        rss <- sum(residuals(rs.lm)^2)
        y <- mf[[attr(terms(lm.object), "response")]]
        syy <- sum(y^2)

        rval <- list(coef = coef(rs.lm), rss = rss, rsquare = (syy - rss) / syy,
                     rstderr = sqrt(rss / (nobs - rs.lm$rank)))

        if(!is.null(new.xpts))
            rval$fitted <- predict(rs.lm, newdata = new.xpts)

        boot.list[[i]] <- rval
    }
    boot.list
}

## R^2 = (SYY - RSS) / SYY

plot.loess.simpleboot <- function(x, add = FALSE, all.lines = FALSE, ...) {
    xpts <- x$new.xpts    
    bmat <- sapply(x$boot.list, "[[", "fitted")
    std <- apply(bmat, 1, sd, na.rm = TRUE)
    orig.pred <- predict(x$orig.loess, data.frame(xpts))
    
    if(!add) {
        xdata <- x$orig.loess$x
        ydata <- x$orig.loess$y
        plot(xdata, ydata, ...)
    }
    lines(xpts, orig.pred)

    if(!all.lines) {
        lines(xpts, orig.pred + 2*std, lty=3)
        lines(xpts, orig.pred - 2*std, lty=3)
    }
    else {
        blist <- x$boot.list
        
        for(i in 1:(x$R)) 
            lines(xpts, blist[[i]]$fitted)
    }
    invisible()          
}

print.loess.simpleboot <- function(x, ...) {
    cat("BOOTSTRAP OF LOESS (method = ", x$method, ")\n\n", sep = "")
    cat("Original Model Fit\n")
    cat("------------------\n")
    print(x$orig.loess)
}

loess.boot <- function(lo.object, R, rows = TRUE, new.xpts = NULL, ngrid = 100, weights = NULL) {
    boot.result <- list()

    if(is.null(new.xpts))
        new.xpts <- seq(min(lo.object$x), max(lo.object$x), len = ngrid)
    if(is.null(weights))
        weights <- rep(1, length(lo.object$x))
    boot.list <- lo.boot.resample(lo.object, R, rows, new.xpts, weights)
    
    boot.result$method <- ifelse(rows, "rows", "residuals")
    boot.result$boot.list <- boot.list
    boot.result$orig.loess <- lo.object
    boot.result$new.xpts <- new.xpts
    boot.result$R <- R
    class(boot.result) <- c("loess.simpleboot")
    boot.result
}

fitted.loess.simpleboot <- function(object, ...) {
    samples(object, "fitted")
}

lo.boot.resample <- function(lo.object, R, rows, new.xpts, weights) {
    f <- formula(lo.object)
    orig.data <- data.frame(lo.object$y, lo.object$x)
    names(orig.data) <- all.vars(f)
    ## Assumes that response is first name and predictor is second name
    
    boot.list <- vector("list", length = R)
    names(boot.list) <- as.character(1:R)

    for(i in 1:R) {
        if(rows) {
            boot.idx <- sample(1:nrow(orig.data), rep = TRUE, prob = weights)
            mf <- orig.data[boot.idx, ]
        }
        else {
            xorig <- lo.object$x
            yorig <- lo.object$y
            res <- yorig - fitted(lo.object)
            boot.res <- sample(res, rep = TRUE, prob = weights)
            yboot <- xorig + boot.res
            mf <- data.frame(yboot, xorig)
            names(mf) <- all.vars(f)
        }
        rs.lo <- loess(f, data = mf, span = lo.object$pars$span)

        rval <- list(rss = sum(residuals(rs.lo)^2),
                     fitted = predict(rs.lo, data.frame(new.xpts)))
        boot.list[[i]] <- rval
    }
    boot.list
}

## lo.boot.res <- function(lo.object, R, new.xpts) {
##     xorig <- lo.object$x
##     yorig <- lo.object$y
##     res <- yorig - fitted(lo.object)
##     f <- formula(lo.object)
## 
##     boot.list <- vector("list", length = R)
##     names(boot.list) <- as.character(1:R)
## 
##     for(i in 1:R) {
##         boot.res <- sample(res, rep = TRUE)
##         yboot <- xorig + boot.res
##         mf <- data.frame(yboot, xorig)
##         names(mf) <- all.vars(f)
##         
##         rs.lo <- loess(f, data = mf, span = lo.object$pars$span)
## 
##         rval <- list(rss = sum(residuals(rs.lo)^2),
##                      fitted = predict(rs.lo, data.frame(new.xpts)))
##         boot.list[[i]] <- rval
##     }
##     boot.list
## }

## lm.boot.rows <- function(lm.object, R, new.xpts) {
##     boot.list <- vector("list", length = R)
##     names(boot.list) <- as.character(1:R)
##     mask.response <- -attr(terms(lm.object), "response")
##     orig.data <- model.frame(lm.object)
##     orig.x <- orig.data[mask.response]
## 
##     for(i in 1:R) {
##         mframe <- model.frame(lm.object)        
##         f <- formula(mframe)
##         boot.idx <- sample(1:nrow(mframe), rep = TRUE)
##         mf <- mframe[boot.idx, ]
##         rs.lm <- lm(f, data = mf)        
## 
##         rval <- list(coef = coef(rs.lm), rss = sum(residuals(rs.lm)^2))
## 
##         rss <- sum(residuals(rs.lm)^2)
##         y <- mf[[attr(terms(rs.lm), "response")]]
##         syy <- sum(y^2)
##         rval$rsquare <- (syy - rss) / syy
##         rval$rstderr <- sqrt(rss / (nrow(mframe) - rs.lm$rank))
##         
##         rval$fitted <- predict(rs.lm, newdata = new.xpts)
## 
##         boot.list[[i]] <- rval
##     }
##     boot.list
## }

## two.boot <- function(sample1, sample2, FUN, R, student = FALSE, M, weights = NULL, ...) {
##     func.name <- ifelse(is.character(FUN), FUN, deparse(substitute(FUN)))
##     func <- match.fun(FUN)
##     ind <- c(rep(1, length(sample1)), rep(2, length(sample2)))
##     data <- list(s1 = sample1, s2 = sample2)
##     extra <- list(...)
##     
##     if(func.name == "quantile") {
##         if(is.na(match("probs", names(extra))))
##             stop("The 'probs' argument must be specified in order to use 'quantile' function")
##         if(length(extra$probs) > 1)
##             stop("Can only bootstrap a single quantile")
##     }
##     boot.func <- function(x) {
##         fval <- func(x$s1, ...) - func(x$s2, ...)
## 
##         if(student) {
##             b <- two.boot(x$s1, x$s2, FUN, R = M, student = FALSE, M = NULL, ...)
##             fval <- c(fval, var(b$t))
##         }
##         fval
##     }
##     ran.gen <- function(d, nothing) {
##         list(s1 = sample(d$s1, rep = TRUE), s2 = sample(d$s2, rep = TRUE))
##     }
##     b <- boot(data, boot.func, R, sim = "parametric", ran.gen = ran.gen)
##     b$data.orig <- b$data
##     b$data <- unlist(data)
##     b$strata <- ind
##     b$sim <- "ordinary"
##     b$stype <- "i"
##     b$student <- student
##     class(b) <- c("simpleboot", class(b))
##     b
## }
