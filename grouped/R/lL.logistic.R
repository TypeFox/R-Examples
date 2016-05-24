"lL.logistic" <-
function(params){
    sigma <- params[p]
    mu <- c(X %*% params[-p])
    -sum(log(plogis(qb, mu, sigma) - plogis(qa, mu, sigma)))
}

