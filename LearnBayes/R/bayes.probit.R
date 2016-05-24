bayes.probit=function (y, X, m, prior = list(beta = 0, P = 0)) 
{
    rtruncated = function(n, lo, hi, pf, qf, ...) qf(pf(lo, ...) + 
        runif(n) * (pf(hi, ...) - pf(lo, ...)), ...)
    if (sum(prior$P)==0) log.marg=NULL
    beta0 = prior$beta
    BI = prior$P
    N = length(y)
    fit = glm(y ~ X - 1, family = binomial(link = probit))
    beta.s = fit$coef
    p = length(beta.s)
    beta = array(beta.s, c(p, 1))
    beta0 = array(beta0, c(p, 1))
    BI = array(BI, c(p, p))
    Mb = array(0, dim = c(m, p))
    lo = c(-Inf, 0)
    hi = c(0, Inf)
    LO = lo[y + 1]
    HI = hi[y + 1]
    postvar=solve(BI + t(X) %*% X)
    aa = chol(postvar)
    BIbeta0 = BI %*% beta0
    post.ord=0
    for (i in 1:m) {
        z = rtruncated(N, LO, HI, pnorm, qnorm, X %*% beta, 1)
        mn = solve(BI + t(X) %*% X, BIbeta0 + t(X) %*% z)
        beta = t(aa) %*% array(rnorm(p), c(p, 1)) + mn
        post.ord=post.ord+dmnorm(beta.s,mn,postvar)
        Mb[i, ] = beta
    }
    if (sum(BI)>0)
    {
    log.f=sum(y*log(fit$fitted)+(1-y)*log(1-fit$fitted))
    log.g=dmnorm(beta.s,beta0,solve(BI),log=TRUE)
    log.marg=log.f+log.g-log(post.ord/m)
    }
    return(list(beta=Mb,log.marg=log.marg))
}
