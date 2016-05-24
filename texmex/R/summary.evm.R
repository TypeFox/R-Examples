summary.evmOpt <- function(object, nsim = 1000, alpha = .050, ...){
    if (ncol(object$data$D$phi) == 1 && ncol(object$data$D$xi) == 1){
        env <- unclass(qqevm(object, nsim = nsim, alpha = alpha))
        env <- list(data=sort(env$dat), envelope=env$sim, Q=env$p)
    }
    else {
        x <- object
        x$data$y <- object$residuals # Standard exponential if GPD model true
        x$threshold <- 0
        x$coefficients <- c(0, 0) # phi not sigma, so 0 not 1

        env <- unclass(qqevm(x, nsim=nsim, alpha=alpha))
        env <- list(data=sort(env$dat), envelope=env$sim, Q=env$p)
    }

    co <- cbind(object$coefficients, object$se, object$coefficients / object$se)
    dimnames(co) <- list(names(coef(object)), c("Value", "SE", "t"))

    res <- list(model = object, coefficients=co, envelope = env,
                nsim = nsim, alpha = alpha)

    oldClass(res) <- "summary.evmOpt"
    res
}

print.summary.evmOpt <- function(x, digits = 3 , ...){
    co <- coef(x)
    env <- x$envelope
    nsim <- x$nsim
    alpha <- x$alpha

    x <- x$model

    cat("Call: ")
    print(x$call, ... )

    cat("\n")
    print(x$family, verbose=TRUE, ...)
    if (is.null(x$penalty) | x$penalty=="none"){
        cat("\nModel fit by maximum likelihood.\n")
    }
    else {
        cat("\nModel fit by penalized maximum likelihood.\n")
    }
    if (x$conv == 0) conv <- TRUE
    else conv <- FALSE
    cat( "\nConvergence:\t\t")
    cat(conv)

    if (x$rate < 1){
        cat( "\nThreshold:\t\t")
        cat(format(unname(x$threshold), digits=digits, ...))
        cat( "\nRate of excess:\t\t")
        cat(format(x$rate, digits=digits, ...))
    }

    cat("\n\nLog-lik.\t\tAIC\n")
    cat(format(x$loglik, digits, ...), "\t\t", format(AIC(x), digits=digits, ...))

    cat( "\n\nCoefficients:\n" )
    print.default(format(co, digits=digits, ...), print.gap=2, quote=FALSE)
    cat( "\n" )

    cat(paste(nsim, " simulated data sets compared against observed data QQ-plot.\n", sep = "" ))
    cat( "Quantile of the observed MSE: ", signif( env$Q, ... ), "\n")
    out <- sum(env$data < env$envelope[1, ] | env$data > env$envelope[2, ])
    level <- paste(round( 100 - alpha*100 , ...), "%", sep = "")
    perc <- round(out / length( env$data ) * 100, digits=digits, ...)
    perc <- paste(perc, "%", sep = "" )
    cat( paste( out, " observations (", perc, ") outside the ", level, " simulated envelope.\n" , sep=""))
    invisible()
}