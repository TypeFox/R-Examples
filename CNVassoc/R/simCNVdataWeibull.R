simCNVdataWeibull<-function (n, mu.surrog, sd.surrog, w, lambda, shape, time.cens=Inf, cnv.random = FALSE)
{
    scale<-lambda^(-1/shape)
    w <- w/sum(w)
    k <- length(mu.surrog)
    if (length(shape)==1)
      shape <- rep(shape, k)
    if (length(scale)==1)
      scale <- rep(scale, k)
    if (cnv.random) {
        nj <- rmultinom(n = 1, size = n, prob = w)
    }
    else {
        nj <- round(n * w)
        nj[1] <- nj[1] + n - sum(nj)
    }
    surrog <- unlist(sapply(seq(along = nj), function(j) rnorm(nj[j],
        mean = mu.surrog[j], sd = sd.surrog[j]), simplify = FALSE))
    cnv <- rep(1:k, nj)
    resp <- unlist(sapply(1:k, function(j) rweibull(nj[j], shape[j], scale[j]), simplify = FALSE))
    cens <- ifelse(resp>time.cens,0,1)
    resp <- ifelse(resp>time.cens,time.cens,resp)    
    out <- data.frame(resp, cens, cnv, surrog)
    out <- out[sample(1:n), ]
    out
}
