copreg=function (x, y, R, S = R, family = 1, exposure = rep(1, length(y)),
    sd.error = FALSE, joint = TRUE, zt = TRUE)
{
    mar <- mle_marginal(x, y, R, S, family, exposure = exposure,
        sd.error = sd.error, zt = zt)
    if (joint == TRUE) {
        joi <- mle_joint(mar$alpha, mar$beta, mar$theta.ifm,
            mar$delta, x, y, R, S, family, exposure, sd.error,
            zt = zt)
    }
    else {
        joi <- mar
    }
    theta_IFM <- mar$theta_IFM
    tau_IFM <- mar$tau_IFM
    alpha <- joi$alpha
    beta <- joi$beta
    delta <- joi$delta
    theta <- joi$theta
    tau <- joi$tau
    sd.alpha <- joi$sd.alpha
    sd.beta <- joi$sd.beta
    family <- joi$family
    ll<-joi$ll
    loglik<-sum(ll)
    alpha0 <- beta0 <- delta0 <- theta0 <- tau0 <- family0 <- sd.alpha0 <- sd.beta0 <- family0 <- loglik0 <-npar0<-
    NULL
    if (joint == TRUE) {
        alpha0 <- mar$alpha
        beta0 <- mar$beta
        delta0 <- mar$delta
        theta0 <- mar$theta
        tau0 <- mar$tau
        sd.alpha0 <- mar$sd.alpha
        sd.beta0 <- mar$sd.beta
        family0 <- mar$family
        loglik0 <- sum(mar$ll)
    }
    npar<-length(alpha)+length(beta)+1
    if (joint==TRUE){
        npar0<-npar
        npar<-npar+1
    }
    outlist <- list(alpha0 = alpha0, beta0 = beta0, delta0 = delta0,
        theta0 = theta0, tau0 = tau0, theta_IFM = theta_IFM,
        tau_IFM = tau_IFM, sd.alpha0 = sd.alpha0, sd.beta0 = sd.beta0,
        family0 = family0, loglik0 = loglik0, npar0=npar0,alpha = alpha,
        beta = beta, delta = delta, theta = theta, tau = tau,
        sd.alpha = sd.alpha, sd.beta = sd.beta, family = family,
        loglik = loglik,npar=npar, zt = zt,ll=ll)
    class(outlist) = "copreg"
    return(outlist)
}
