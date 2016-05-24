gpdRangeFit <-
function (data, umin=quantile(data, .05), umax=quantile(data, .95),
          nint = 10, 
          penalty="gaussian", priorParameters=NULL, alpha=.05) {

    m <- s <- hi <- lo <- matrix(0, nrow = nint, ncol = 2)
    u <- seq(umin, umax, length = nint)
    qz <- qnorm(1-alpha/2)
    for (i in 1:nint) {
        z <- evm(data, th=u[i], penalty=penalty, priorParameters=priorParameters)
        m[i, ] <- z$coefficients
        m[i, 1] <- m[i, 1] - m[i, 2] * u[i]
        d <- matrix(c(1, -u[i]), ncol = 1)
        v <- t(d) %*% z$cov %*% d
        s[i, ] <- sqrt(diag(z$cov))
        s[i, 1] <- sqrt(v)
        
        hi[i, ] <- m[i, ] + qz * s[i, ]
        lo[i, ] <- m[i, ] - qz * s[i, ]
    }
    res <- list(th=u, par=m , hi=hi, lo=lo, data=data)
    oldClass(res) <- 'gpdRangeFit'
    res
}

print.gpdRangeFit <- function(x, ...){
    sc <- cbind(threshold=x$th, phi=x$par[, 1], lo=x$lo[, 1], hi=x$hi[, 1])
    sh <- cbind(threshold=x$th, xi=x$par[, 2], lo=x$lo[, 2], hi=x$hi[, 2])

    print(sc); print(sh)
    invisible()
}

summary.gpdRangeFit <- function(object, ...){
    sc <- cbind(threshold=object$th, phi=object$par[, 1], lo=object$lo[, 1], hi=object$hi[, 1])
    sh <- cbind(threshold=object$th, xi=object$par[, 2], lo=object$lo[, 2], hi=object$hi[, 2])

    list(phi=summary(sc), xi=summary(sh))
}

plot.gpdRangeFit <- function(x, xlab="Threshold", ylab=NULL,
                             main=NULL, addNexcesses=TRUE, ...){
    #############################################################
    ## Get axis labels and titles
    if (missing(ylab)){
        ylab <- c("log(scale)", "shape")
    }
    else if (length( ylab ) != 2){
        stop("length of ylab should be 2")
    }

    if (!missing(main) && length(main) != 2){
        stop("length of main should be 2")
    }
    ##
    #############################################################

    data <- x$data

    for (i in 1:2) {
        yl <- range(x$hi[, i], x$lo[, i])
        plot(x$th, x$par[, i], ylim = yl, type = "b",
             xlab=xlab, ylab=ylab[i], main=main[i], ...)
        for (j in 1:length(x$th)){
            lines(c(x$th[j], x$th[j]), c(x$hi[j, i], x$lo[j, i]))
        }
        if(addNexcesses){
          axis(3,at=axTicks(1),labels=sapply(axTicks(1),function(u)sum(data>u)),cex=0.5)
          mtext("# threshold excesses")
        }
    }
    invisible()
}

