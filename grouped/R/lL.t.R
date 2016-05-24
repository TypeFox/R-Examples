"lL.t" <-
function(params){
    sigma <- params[p]
    mu <- c(X %*% params[-p])
    -sum(log(pt((qb - mu) / sigma, df.) - pt((qa - mu) / sigma, df.)))
}

