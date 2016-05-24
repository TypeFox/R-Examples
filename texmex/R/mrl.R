mrl <-
function (data, umin = min(data), umax = max(data) - 0.1,
          nint = 100, alpha=.050) {
    data <- c(data)
    AllData <- data
    x <- lo <- hi <- numeric(nint)
    Threshold <- seq(umin, umax, length = nint)
    z <- qnorm(1 - alpha/2)
    for (i in 1:nint) {
        data <- data[data > Threshold[i]]
        x[i] <- mean(data - Threshold[i])
        sdev <- sqrt(var(data))
        lo[i] <- x[i] - z*sdev/sqrt(length(data))
        hi[i] <- x[i] + z*sdev/sqrt(length(data))
    }
    
    res <- cbind(threshold=Threshold, MRL=x, lo=lo, hi=hi)
    res <- list(mrl=res, data=AllData)
    oldClass(res) <- 'mrl'
    res
}

print.mrl <- function(x, ...){
    print(x$mrl)
    invisible()
}

summary.mrl <- function(object, ...){
    summary(object$mrl)
}

plot.mrl <- function(x, xlab="Threshold", ylab="Mean excess", ...){

    data <- x$data
    x <- x$mrl

    th <- x[, "threshold"]
    mrl <- x[, "MRL"]
    xl <- x[, "lo"]
    xu <- x[, "hi"]

    plot(th, mrl, type = "l", xlab = xlab, ylab = ylab,
        ylim = c(min(xl[!is.na(xl)]), max(xu[!is.na(xu)])), ...)
    lines(th[!is.na(xl)], xl[!is.na(xl)], lty = 2)
    lines(th[!is.na(xu)], xu[!is.na(xu)], lty = 2)
    rug(data)
    invisible()
}

