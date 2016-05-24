bayes.model.selection=function (y, X, c, constant = TRUE) 
{
    base2 = function(s, k) {
        r = rep(0, k)
        for (j in seq(k, 1, by = -1)) {
            r[j] = floor(s/(2^(j - 1)))
            s = s - r[j] * (2^(j - 1))
        }
        return(r)
    }
    regpost.mod = function(theta, stuff) {
        y = stuff$y
        X = stuff$X
        c = stuff$c
        beta = theta[-length(theta)]
        sigma = exp(theta[length(theta)])
        if (length(beta) > 1) 
            loglike = sum(dnorm(y, mean = X %*% as.vector(beta), 
                sd = sigma, log = TRUE))
        else loglike = sum(dnorm(y, mean = X * beta, sd = sigma, 
            log = TRUE))
        logprior = dmnorm(beta, mean = 0 * beta, varcov = c * 
            sigma^2 * solve(t(X) %*% X), log = TRUE)
        return(loglike + logprior)
    }
    require(LearnBayes)
    X = as.matrix(X)
    if (constant == FALSE) 
        X = cbind(1, X)
    p = dim(X)[2] - 1
    GAM = array(TRUE, c(2^p, p + 1))
    for (k in 1:(2^p)) GAM[k, ] = as.logical(c(1, base2(k - 1, 
        p)))
    gof = rep(0, 2^p)
    converge = rep(TRUE, 2^p)
    for (j in 1:2^p) {
        X0 = X[, GAM[j, ]]
        fit = lm(y ~ 0 + X0)
        beta = fit$coef
        s = sqrt(sum(fit$residuals^2)/fit$df.residual)
        theta = c(beta, log(s))
        S = list(X = X0, y = y, c = c)
        fit = laplace(regpost.mod, theta, S)
        gof[j] = fit$int
        converge[j] = fit$converge
    }
    Prob=exp(gof-max(gof))/sum(exp(gof-max(gof)))
    mod.prob=data.frame(GAM[, -1], round(gof,2), round(Prob,5))
    names(mod.prob)=c(dimnames(X)[[2]][-1],"log.m","Prob")
    return(list(mod.prob=mod.prob, converge = converge))
}
