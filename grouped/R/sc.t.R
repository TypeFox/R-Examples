"sc.t" <-
function(params){
    sigma <- params[p]
    mu <- c(X %*% params[-p])
    fnb <- dgt(qb, mu, sigma, df.)
    fna <- dgt(qa, mu, sigma, df.)
    dF <- pt((qb - mu) / sigma, df.) - pt((qa - mu) / sigma, df.)
    sbeta <- colSums( X * (fnb - fna) / dF )
    ssigma <- sum( (fnb * (qb. - mu) - fna * (qa. - mu)) / (sigma * dF) )
    c(sbeta, ssigma)
}

