evmBoot <- function(o, R=1000, trace=100, theCall){
    if (class(o) != "evmOpt"){
        stop("o must be of class 'evmOpt'")
    }

    if (missing(theCall)){ theCall <- match.call() }

    d <- o$data
    param <- texmexMakeParams(coef(o), d$D)
    rng <- o$family$rng

    bfun <- function(i){
        if (i %% trace == 0){ cat("Replicate", i, "\n") }

        d$y <- rng(length(d$y), param, o)

        evmFit(d, o$family, th=o$threshold, prior=o$penalty,
               priorParameters=o$priorParameters,
               start=o$coefficients, hessian=FALSE)$par
    }

    res <- t(sapply(1:R, bfun))

    se <- apply(res, 2, sd)
    b <- apply(res, 2, mean) - coef(o)

    if (any(abs(b/se) > .25)){
        warning("Ratio of bias to standard error is high")
    }

    res <- list(call=theCall, replicates=res, map=o)
    oldClass(res) <- "evmBoot"
    res
}

print.evmBoot <- function(x, ...){
    print(x$call)
    means <- apply(x$replicates, 2, mean)
    medians <- apply(x$replicates, 2, median)
    sds <- apply(x$replicates, 2, sd)
    bias <- means - x$map$coefficients
    res <- rbind(x$map$coefficients, means, bias, sds, medians)
    rownames(res) <- c("Original", "Bootstrap mean", "Bias", "SD", "Bootstrap median")
    print(res, ...)
    if (any(abs(res[3,] / res[4,]) > .25)){
        warning("Ratio of bias to standard error is high")
    }
    invisible(res)
}

coefficients.evmBoot <- coef.evmBoot <- function(object, ...){
    apply(object$replicates, 2, mean)
}

summary.evmBoot <- function(object, ...){
    means <- apply(object$replicates, 2, mean)
    medians <- apply(object$replicates, 2, median)
    sds <- apply(object$replicates, 2, sd)
    bias <- means - coef(object$map)
    res <- rbind(coef(object$map), means, bias, sds, medians)
    rownames(res) <- c("Original", "Bootstrap mean", "Bias", "SD", "Bootstrap median")

    if (any(abs(res[3,] / res[4,]) > .25)){
        warning("Ratio of bias to standard error is high")
    }

covs <- var(object$replicates)
    res <- list(call = object$call, margins=res, covariance=covs)
    oldClass(res) <- "summary.evm.boot"
    res
}

print.summary.evmBoot <- function(x, ...){
    print(x$call)
    print(x$margins)
    cat("\nCorrelation:\n")
    print(cov2cor(x$covariance))
    invisible()
}

plot.evmBoot <- function(x, col=4, border=FALSE, ...){
    pfun <- function(x, col, border, xlab,...){
        d <- density(x, n=100)
        hist(x, prob=TRUE, col=col, border=border, main="", xlab=xlab, ...)
        lines(d, lwd=2, col="grey")
        rug(x)
        invisible()
    }
    for (i in 1:ncol(x$replicates)){
        pfun(x$replicates[,i], xlab=colnames(x$replicates)[i], col=col, border=border)
        abline(v=coef(x$map)[i], col="cyan")
    }
    invisible()
}

