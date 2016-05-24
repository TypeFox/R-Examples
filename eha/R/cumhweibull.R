Hweibull <- function(x, shape, scale = 1, log = FALSE){
    if (shape <= 0 || scale <= 0)
      stop("scale and shape must be positive")
    res <- ifelse(x < 0, 0, (x / scale)^shape) 
    if (log) res <- logb(res, base = exp(1))
    
    res
}
