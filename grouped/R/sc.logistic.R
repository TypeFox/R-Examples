"sc.logistic" <-
function(params){
    sigma <- params[p]
    mu <- c(X %*% params[-p])
    fnb <- dlogis(qb, mu, sigma)
    fna <- dlogis(qa, mu, sigma)
    dF <- plogis(qb, mu, sigma) - plogis(qa, mu, sigma)
    sbeta <- colSums( X * (fnb - fna) / dF )
    ssigma <- sum( (fnb * (qb. - mu) - fna * (qa. - mu)) / (sigma * dF) )
    c(sbeta, ssigma)
}

