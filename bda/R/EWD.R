pewd <- function(q, alpha, shape, scale){
    (pweibull(q, shape=shape, scale=scale))^alpha
}

dewd <- function(x, alpha, shape, scale){
    fx <- dweibull(x, shape=shape, scale=scale)
    Fx <- pweibull(x, shape=shape, scale=scale)
    alpha * Fx^(alpha-1) * fx 
}

qewd <- function(p, alpha, shape, scale){
    (-log(1-p^(1/alpha)))^(1/shape)*scale
}
