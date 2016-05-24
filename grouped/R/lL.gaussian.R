"lL.gaussian" <-
function(params){
    sigma <- params[p]
    mu <- c(X %*% params[-p])
    -sum(log(pnorm(qb, mu, sigma) - pnorm(qa, mu, sigma)))
}

