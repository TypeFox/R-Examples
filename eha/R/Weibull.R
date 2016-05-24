hweibull <- function(x, shape, scale = 1, log = FALSE){
    if (any(shape <= 0) || any(scale <= 0))
      stop("scale and shape must be positive")
    res <- ifelse(x < 0, 0, shape * (x / scale)^(shape - 1) / scale)
    if (log) res <- log(res)
    res
}

Hweibull <- function(x, shape, scale = 1, log.p = FALSE){
    if (any(shape <= 0) || any(scale <= 0))
      stop("scale and shape must be positive")
    res <- ifelse(x < 0, 0, (x / scale)^shape)
    if (log) res <- logb(res, base = exp(1))

    res
}

