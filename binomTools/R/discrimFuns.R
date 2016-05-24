##################################################################
### Rune Haubo Bojesen Christensen, rhbc@imm.dtu.dk
### August 5th, 2011
###
### This script contains R-functions to estimate the performance
### measures described in Tjur T. (2009). Coefficients of determination in
### logistic regression models - a new proposal: The coefficient of
### discrimination. The American Statistician, vol 63, no. 4: 366-372
### (http://staff.cbs.dk/tuetjur/R2.pdf)
###
### As an example copy the following commands to the R-promt:
###     ## Get example of a binomial glm:
###     example(predict.glm)
###
###     ## Obtain access to discrim-functions:
###     source("discrimFuns.R")
###
###     ## Test Rsq function and its print and plot methods:
###     (Rsq.budworm <- Rsq(budworm.lg))
###     HLtest(Rsq.budworm)
###     HLtest(Rsq.budworm, method="fixed")
###     X2GOFtest(Rsq.budworm)
###     plot(Rsq.budworm, "hist") ## or simply 'plot(Rsq.budworm)'
###     plot(Rsq.budworm, "ecdf")
###     plot(Rsq.budworm, "ROC")
###
###
###  This program is free software; you can redistribute it and/or modify
###  it under the terms of the GNU General Public License as published by
###  the Free Software Foundation; either version 2 or 3 of the License
###  (at your option).
###
###  This program is distributed in the hope that it will be useful,
###  but WITHOUT ANY WARRANTY; without even the implied warranty of
###  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###  GNU General Public License for more details.
###
###  A copy of the GNU General Public License is available at
###  http://www.r-project.org/Licenses/
###
##################################################################

Rsq <- function(object, ...) {
  UseMethod("Rsq")
}

##eAUC <- function(object, ...) {
##  UseMethod("eAUC")
##}

HLtest <- function(object, ...) {
  UseMethod("HLtest")
}

X2GOFtest <- function(x, ...) {
    UseMethod("X2GOFtest")
}

group <- function(object, ...) {
    UseMethod("group")
}

Rsq.glm <- function(object, ...) {
### Initial testing:
    stopifnot(object$family$family == "binomial")
###    stopifnot(object$family$link == "logit")
### It might be useful even with other links than the logit although D
### will have another interpretation.
    y <- object$y
    N <- length(y)
    wt <- object$prior.weights
    muHat <- object$fitted.values
    resp <- cbind(suc = y*wt , n = wt)
    Y <- apply(resp, 1,
               function(x) c(rep(1, x[1]), rep(0, x[2]-x[1])))
    if(is.list(Y))
        Y <- unlist(Y)
    else
        Y <- c(Y)
    piHat <- c(apply(cbind(muHat, wt), 1,
                     function(x) rep(x[1], x[2])))
    if(is.list(piHat))
        piHat <- unlist(piHat)
    else
        piHat <- c(piHat)

### Sums of squares:
    SSDres <- sum((Y-piHat)^2)
    SSDmod <- sum((piHat - mean(Y))^2)
    SSDtot <- sum((Y - mean(Y))^2)

### R-square measures and the coefficient of discrimination:
    R2mod <- SSDmod/SSDtot
    R2res <- 1 - SSDres/SSDtot
    R2cor <- cor(Y, piHat)^2
    D <- (R2mod + R2res)/2
    ## D <- sqrt(R2mod*R2cor) is another possibility

    res <- list(Y=Y, piHat=piHat, N=N, glm=object,
                R2mod=R2mod, R2res=R2res, R2cor=R2cor, D=D)

    class(res) <- "Rsq"
    res
}

print.Rsq <- function(x, digits = getOption("digits"), ...) {
    res <- with(x, data.frame(R2mod=R2mod, R2res=R2res, R2cor=R2cor, D=D))
    cat("\nR-square measures and the coefficient of discrimination, 'D':\n\n")
    print.data.frame(res, row.names=FALSE, digits=digits)
    cat("\nNumber of binomial observations: ", x$N)
    cat("\nNumber of binary observation: ", length(x$Y))
    cat("\nAverage group size: ", with(x, length(Y)/N), "\n\n")
    invisible(x)
}

plot.Rsq <- function(x, which=c("hist", "ecdf", "ROC"), ...) {
    which <- match.arg(which)
    pi0 <- with(x, piHat[Y == 0])
    pi1 <- with(x, piHat[Y == 1])
    F0 <- ecdf(pi0)
    F1 <- ecdf(pi1)
    oldPar <- par(no.readonly=TRUE)
    on.exit(par(oldPar))

    if(which == "hist") {
        breaks <- seq(0, 1, .1)
        par(mfrow=c(2,1), mar=c(4,4,0,0)+.1)
        x$h0 <- hist(x$piHat[x$Y == 0], breaks=breaks, main="",
                     xlab=expression(hat(pi)['y=0']))
        yrange <- range(x$h0$counts)
        x$h1 <- hist(x$piHat[x$Y == 1], breaks=breaks, ylim=yrange, main="",
                     xlab=expression(hat(pi)['y=1']))
    }
    if(which == "ecdf") {
        par(mfrow=c(2,1), mar=c(3,3,0,0)+.5)
        piVals <- seq(0, 1, length=1000)
        plot(piVals, F0(piVals), type="l", xlab="", ylab="")
        mtext(expression(pi), side=1, line=2)
        mtext(expression(F['y=0'](pi)), side=2, line=2)
        plot(piVals, F1(piVals), type="l", xlab="", ylab="")
        mtext(expression(pi), side=1, line=2)
        mtext(expression(F['y=1'](pi)), side=2, line=2)
    }
    if(which == "ROC") {
        par(mfrow=c(1,1), mar=c(4,4,0,0)+.5)
        piVals <- seq(0, 1, length=1000)
        plot(1 - F0(piVals), 1 - F1(piVals), type="l",
             xlab=expression(1 - F['y=0'](pi)),
             ylab=expression(1 - F['y=1'](pi)))
        lines(c(0,1), c(0,1), lty=2)
    }

    class(x) <- "plot.Rsq"
    invisible(x)
}

##eAUC.Rsq <- function(x, ...) {
##    pi0 <- with(x, piHat[Y == 0])
##   pi1 <- with(x, piHat[Y == 1])
##    F0 <- ecdf(pi0)
##    F1 <- ecdf(pi1)
##    piVals <- seq(0, 1, length=1000)
##    F0inv <- approxfun(piVals, F0(1 - piVals), method = "constant")
##    intFn <- function(x) 1 - F1(F0inv(1 - x))
##    integrate(f = intFn, lower = 0, upper = 1)
##}

HLtest.Rsq <- function(object, method=c("deciles", "fixed"), decile.type=8,
                       ...) {
    method <- match.arg(method)
    piHat <- object$piHat
    Y <- object$Y
    ## Indicator variables of group attribution
    if(method == "deciles") {
        q0 <- quantile(piHat, probs=seq(.1, .9, 0.1),
                       type=decile.type)
        i0 <- cut(piHat[Y == 0], breaks=c(0, q0, 1), labels=1:10)
        i1 <- cut(piHat[Y == 1], breaks=c(0, q0, 1), labels=1:10)
        i <- cut(piHat, breaks=c(0, q0, 1), labels=1:10)
    }
    if(method == "fixed") {
        i0 <- cut(piHat[Y == 0], breaks=seq(0, 1, .1), labels=1:10)
        i1 <- cut(piHat[Y == 1], breaks=seq(0, 1, .1), labels=1:10)
        i <- cut(piHat, breaks=seq(0, 1, .1), labels=1:10)
    }
    ## The sum of the probabilities in each of the 2x10 groups,
    ## i.e. the expected number of observations in each group:
    E1 <- tapply(piHat, i, sum)
    E <- rbind(table(i)-E1, E1)
    ## Observed no. observations in each group
    O <- matrix(0, 2, 10)
    ti0 <- table(i0)
    ti1 <- table(i1)
    O[1,as.numeric(names(ti0))] <- ti0
    O[2,as.numeric(names(ti1))] <- ti1
    colnames(O) <- colnames(E)

    resid <- (O-E)/sqrt(E)
    E[is.na(E)] <- 0
    resid[is.na(resid)] <- 0
    X2 <- sum(resid^2, na.rm=TRUE)
    p.value <- pchisq(X2, 10-2, lower=FALSE)
    res <- list(expected=E, observed=O, resid=resid, X2=X2,
                p.value=p.value, method=method)
    class(res) <- "HLtest.Rsq"
    res
}

print.HLtest.Rsq <- function(x, digits = getOption("digits"), ...) {
    cat("\nHosmer and Lemeshow's goodness-of-fit test for logistic regression models\n")
    cat("with method =", x$method)
    cat("\nObserved counts:\n")
    print.data.frame(as.data.frame(x$observed), row.names=FALSE)
    cat("\nExpected counts:\n")
    print.data.frame(as.data.frame(round(x$expected,1)), row.names=FALSE)
    cat("\nPearson residuals:\n")
    print.data.frame(as.data.frame(round(x$resid,1)), row.names=FALSE)
    cat("\nChi-square statistic: ", x$X2, " with ", 10-2, " df\n")
    cat("P-value: ", x$p.value, "\n\n")

    if(any(x$expected < 5))
        cat("Warning: Test may not be appropriate due to low expected frequencies\n\n")
    invisible(x)
}

X2GOFtest.Rsq <- function(x, ...) {
    obj <- x$glm
    ##    piHat <- x$piHat
    if(any(duplicated(model.matrix(obj))))
        warning("model contains duplicated covariate patterns, so the test may not be accurate")
    piHatJ <- fitted(obj) ## unique(x$piHat)
    wt <- obj$prior.weight
    yy <- obj$y * wt
    vJ <- wt * piHatJ * (1 - piHatJ)
    ##    names(vJ) <- NULL
    cJ <- (1 - 2 * piHatJ) / vJ
    X2 <- sum(resid(obj, type="pearson")^2)
    ## ind <- match(piHatJ, fitted(obj))

    form <- as.formula(obj$formula)
    form[[2]] <- as.name("cJ") ## replace response variable
    dat <- obj$data
    dat$cJ <- cJ
    dat$vJ <- vJ
    RSS <- sum(resid(lm(form, data=dat, weights=vJ))^2)
    A <- 2*(length(piHatJ) - sum(1/wt))
    z <- (X2 - obj$df.residual)/sqrt(A + RSS)
    p.value <- 2 * pnorm(abs(z), lower=FALSE)

    res <- list(p.value=p.value, z.score=z, RSS=RSS, X2=X2)
    class(res) <- "X2GOFtest.Rsq"
    res
}

print.X2GOFtest.Rsq <- function(x, ...) {
    cat("\nPearson goodness-of-fit test for logistic regression models
with normal approximation\n")
    cat("\nChi-square statistic: ", x$X2, "\n")
    cat("Normal deviate", x$z.score, "\n")
    cat("P-value: ", x$p.value, "\n\n")
    invisible(x)
}

group.glm <- function(object, eval=TRUE, ind=NULL, ...) {
### or should it be aggregate.glm?

### Should the model allow for additional variables to be taken
### account of? e.g. 'extra.var=NULL'

    if(is.null(object$model))
        stop("glm object has to be fitted with 'model=TRUE'")
    if(object$family$family != "binomial")
        stop("Method only implemented for binomial glms")

    oldData <- object$model
    oldN <- nrow(oldData)
    y <- oldData[[1]] # Should probably be:
##    y <- model.response(object$model)

    if(NCOL(y) == 1L) {
        cnam <- c("successes", "failures")
        if (is.factor(y)) {
            s <- y != levels(y)[1L]
            cnam <- levels(y)[1:2]
        }
        if (any(y < 0 | y > 1))
            stop("y values must be 0 <= y <= 1")
        s <- y
        f <- rep.int(1, oldN) - s
    } else if(NCOL(y) == 2L) {
        s <- y[,1]
        f <- y[,2]
        cnam <- colnames(y)
    } else stop("reponse not recognized")

    if(!is.null(ind)){
        ind <- unclass(ind)
        if(oldN != length(ind))
            stop("length(ind) has to match length(object$model)")
        keep <- !duplicated(ind) ## row indicator
    } else { ## construct 'ind' from covariate pattern in data:

### Note that problems may occur with some usage of 'offset',
### i.e. when offset is used both in formula and as separate
### argument.

        varnames <- sapply(attr(object$terms, "variables"),
                           deparse)[-(1:2)]
        X <- object$model[varnames]
         if(!is.null(model.offset(object$model)))
            X['(offset)'] <- model.offset(object$model)
        oldRows <- apply(X, 1, function(x) paste(x, collapse = "\r"))
        keep <- !duplicated(X)
        uniqueX <- X[keep, , drop=FALSE]
        newRows <- apply(uniqueX, 1, function(x) paste(x, collapse = "\r"))
        ind <- match(oldRows, newRows)
        newN <- length(newRows)
    }

### ? Check if supplied 'ind' is a subset of covariate pattern?
### If it is not, the user is warned that original and new models are
### not identical.

    newN <- sum(keep)
    newData <- oldData[keep, -1, drop=FALSE]

    if(!is.null(model.weights(oldData))){ ## Extract original weights
        wt <- model.weights(oldData)
        newData <- newData[-match("(weights)", names(newData))]
        if (any(wt%%1 != 0))
            warning("Aggregation might not be appropriate with non-integer weights")
    } else wt <- rep.int(1, oldN)

    newData$YY <- cbind(tapply(s*wt, ind, sum),
                        tapply(f*wt, ind, sum))
    colnames(newData$YY) <- cnam

    lcall <- as.list(update(object, formula=YY ~ ., data=newData,
                            evaluate=FALSE))
    if("subset" %in% names(lcall)) ## remove 'subset' if it exists
        lcall <- lcall[-match("subset", names(lcall))]
    if("weights" %in% names(lcall)) ## remove 'weights' if it exists
        lcall <- lcall[-match("weights", names(lcall))]
    call <- as.call(lcall)

### ? Should old object be included in output?

    res <- list(newCall=call, newData=newData, oldData=oldData,
                oldN=oldN, newN=newN, oldObject=object)
    if(eval) {
        newObject <- eval(call)
        res$newObject <- newObject
        res$equal <- identical(round(coef(object), 3),
                               round(coef(newObject), 3))
        if(!res$equal)
            warning("Original and new models appear not to be identical")
    }
    class(res) <- "group.glm"
    res
}
