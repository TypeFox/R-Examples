blinreg=function (y, X, m, prior=NULL) 
{
    if(length(prior)>0)
      { c0=prior$c0; beta0=matrix(prior$b0,c(1,length(prior$b)))}

    fit = lm(y ~ 0 + X)
    bhat = matrix(fit$coef, c(1, fit$rank))
    s2 = sum(fit$residuals^2)/fit$df.residual

    if(length(prior)==0)
    {
    shape = fit$df.residual/2
    rate = fit$df.residual/2 * s2
    beta.m = bhat
    vbeta = vcov(fit)/s2
    } else
    {
    shape = length(y)/2
    rate = fit$df.residual/2 * s2 + 
         (beta0 - bhat) %*% t(X) %*% X %*% t(beta0 - bhat)/2/(c0+1)
    beta.m = c0/(c0+1)*(beta0/c0 + bhat)
    vbeta = vcov(fit)/s2*c0/(c0+1)
    }

    sigma = sqrt(1/rgamma(m, shape = shape, rate = rate))
    beta = rmnorm(m, mean=rep(0, fit$rank), varcov=vbeta) 
    beta = array(1, c(m, 1)) %*% beta.m +
           array(sigma, c(m, fit$rank))*beta

    return(list(beta = beta, sigma = sigma))
}

