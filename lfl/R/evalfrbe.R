#.u2 <- function(predicted, real, rw) {
    #sqrt(sum((real - predicted)^2) / length(real)) / sqrt(sum((real - rw)^2) / length(real))
#}


evalfrbe <- function(fit, 
                     real,
                     error=c('smape', 'mase', 'rmse')) {
    if (!is.frbe(fit)) {
        stop("'fit' must be an instance of class 'frbe'")
    }
    if (!is.vector(real) || !is.numeric(real)) {
        stop("'real' must be a numeric vector")
    }
    error <- match.arg(error)
    if (error == 'smape') {
        errorFunc <- smape
    } else if (error == 'mase') {
        errorFunc <- mase
    } else if (error == 'rmse') {
        errorFunc <- rmse
    } else {
        stop("Unknown error function name")
    }

    d <- fit$forecasts
    d$avg <- rowSums(d) / ncol(d)
    d$frbe <- fit$mean

    if (length(real) > nrow(d)) {
        length(real) <- nrow(d)
    } else {
        d <- d[seq_along(real), ]
    }

    r <- colwise(function(col) { errorFunc(col, real) })(d)
    #u2 <- colwise(function(col) { .u2(col, real, fit$forecasts$randomWalk) })(d)
    #names(u2) <- paste('u2', names(u2))
    #return(cbind(r, u2))
    return(r)
}
