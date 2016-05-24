### The Gompertz-Makeham distribution, the hazard is
### h(t; shape, scale) = shape[2] + shape[1] * exp(t / scale), t >= 0

pmakeham <- function(q, shape = c(1, 1), scale = 1,
                     lower.tail = TRUE, log.p = FALSE){

    if (length(shape) != 2) stop("shape must have length 2")
    if (any(c(shape, scale) <= 0)) stop("shape and scale must be positive")

    if (all(q <= 0)){ # Fix trivial case first
        ret <- 1 - as.numeric(lower.tail)
        if (log.p) ret <- log(ret)
        return ( rep(ret, length(q)) )
    }
    ## Serious business:

    tmp <- ifelse(q <= 0,
                  0,
                  -( shape[2] * q + shape[1] * scale * expm1(q / scale) )
                  )
    if (lower.tail){
        if (log.p){
            ret <- ifelse(q <= 0, -Inf, log1p(-exp(tmp)))
        }else{
            ret <- ifelse(q <= 0, 0, -expm1(tmp))
        }
    }else{
        if (log.p){
            ret <- ifelse(q <= 0, 0, tmp)
        }else{
            ret <- ifelse(q <= 0, 1, exp(tmp))
        }
    }

    return ( ret )
}

hmakeham <- function(x, shape = c(1, 1), scale = 1, log = FALSE){

    if (length(shape) != 2) stop("shape must have length 2")
    if (any(c(shape, scale) <= 0)) stop("shape and scale must be positive")

    ret <- ifelse(x < 0,
                  0,
                  shape[2] + shape[1] * exp(x / scale))
    if (log) ret <- log(ret)

    return ( ret )
}

dmakeham <- function(x, shape = c(1, 1), scale = 1, log = FALSE){

    if (length(shape) != 2) stop("shape must have length 2")
    if (any(c(shape, scale) <= 0)) stop("shape and scale must be positive")

    ret <- ifelse(x < 0, 0, hmakeham(x, shape, scale, log = FALSE) *
                  pmakeham(x, shape, scale,
                           lower.tail = FALSE, log.p = FALSE))
    if (log) ret <- log(ret)

    return ( ret )
}

Hmakeham <- function(x, shape = c(1, 1), scale = 1, log.p = FALSE){

    if (length(shape) != 2) stop("shape must have length 2")
    if (any(c(shape, scale) <= 0)) stop("shape and scale must be positive")

    ret <- ifelse(x <= 0, 0,
                  shape[2] * x + shape[1] * scale * expm1(x / scale))
    if (log.p) ret <- log(ret)

    return ( ret )
}

qmakeham <- function(p, shape = c(1, 1), scale = 1,
                     lower.tail = TRUE, log.p = FALSE){

    if (length(shape) != 2) stop("shape must have length 2")
    if (any(c(shape, scale) <= 0)) stop("shape and scale must be positive")

    stop("Sorry, qmakeham is not yet implemented")
}

rmakeham <- function(n, shape = c(1, 1), scale = 1){

    if (length(shape) != 2) stop("shape must have length 2")
    if (any(c(shape, scale) <= 0)) stop("shape and scale must be positive")

    ## In the absence of 'qmakeham', we utilize the fact that a
    ## Gomperts-Makeham distributed random variable has the same
    ## distribution as the minimum of two independent random variables,
    ## one with an exponential and the other with a Gompertz distribution.

    x1 <- rexp(n, rate = shape[2])
    x2 <- rgompertz(n, shape = shape[1], scale = scale)

    return ( pmin(x1, x2) )
}
