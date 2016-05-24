### Contains functions for the Gompertz distribution.

pEV <- function(q, shape = 1, scale =  1,
                      lower.tail = TRUE, log.p = FALSE){
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    q <- ifelse(q <= 0, q, (q / scale)^shape)
    if (log.p){
        if (lower.tail){
            ret <- ifelse(q <= 0, -Inf, log1p(-exp(-expm1(q))))
        }else{
            ret <- ifelse(q <= 0, 0, -expm1(q))
        }
    }else{
        if (lower.tail){
            ret <- ifelse(q <= 0, 0, 1.0 - exp(-expm1(q)))
        }else{
            ret <- ifelse(q <= 0, 1, exp(-expm1(q)))
        }
    }

    return ( ret )
}

dEV <- function(x, shape = 1, scale = 1, log = FALSE){
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    y <- (x / scale)^shape

    ret <- ifelse(x < 0, 0, (shape / scale) * (x / scale)^(shape - 1) *
                  exp(y - expm1(y)))

    if (log) ret <- log(ret)

    return ( ret )
}

hEV <- function(x, shape = 1, scale = 1, log = FALSE){
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    y <- (x / scale)^shape

    ret <- ifelse(x < 0, 0,
                  (shape / scale) * (x / scale)^(shape - 1) * exp(y))

    if (log) ret <- log(ret)

    return ( ret )
}

qEV <- function(p, shape = 1, scale = 1,
                      lower.tail = TRUE, log.p = FALSE){
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p

    ok <- (p >= 0) & (p <= 1)

    ret <- ifelse(ok, (1 / scale) * (log1p(-log1p(-p)))^(1 / shape), NaN)

    if (!all(ok)) warning("qEV produced NaN's")

    return ( ret )
}

HEV <- function(x, shape = 1, scale = 1, log.p = FALSE){
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    y <- ifelse(x <= 0, 0, (x / scale)^shape)

    if (log.p){
        ret <- ifelse(y <= 0, -Inf, log(expm1(y)))
    }else{
        ret <- ifelse(y <= 0, 0, expm1(y))
    }

    return ( ret )
}

rEV <- function(n, shape = 1, scale = 1){
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    y <- runif(n)

    return ( (1 / scale) * (log1p(-log1p(-y)))^(1 / shape) )
}





