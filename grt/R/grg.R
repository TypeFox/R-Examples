grg <- function(response, fixed = FALSE, k = 2)
{
    response <- as.vector(response)
    N <- length(response)
    n <- table(response)
    if(length(n) > 2) stop("more than 2 categories found")
    fit <- NULL
    if(fixed){
        fit$logLik <- N * log(fit$par <- .5)
        attr(fit$logLik, "df") <- 0
    } else {
        m <- max(n)
        fit$par <- p <- m/N
        fit$logLik <- m*log(p) + (N-m)*log(1-p)
        if(p==1) fit$logLik <- m*log(p)
        attr(fit$logLik, "df") <- 1
    }
    class(fit$logLik) <- "logLik"
    fit$AIC <- AIC(fit$logLik, k = k)
    class(fit) <- "grg"
    fit
}

extractAIC.grg <- function(fit, scale, k = 2, ...)
{
    loglik <- fit$logLik
    edf <- attr(loglik, "df")
    c(edf, -2 * loglik + k * edf)
}