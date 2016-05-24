"sc.gaussian" <-
function(params){
    sigma <- params[p]
    mu <- c(X %*% params[-p])
    fnb <- dnorm(qb, mu, sigma)
    fna <- dnorm(qa, mu, sigma)
    dF <- pnorm(qb, mu, sigma) - pnorm(qa, mu, sigma)
    sbeta <- colSums(X * (fnb - fna) / dF)
    ssigma <- sum((fnb * (qb. - mu) - fna * (qa. - mu)) / (sigma * dF))
    c(sbeta, ssigma)
}

