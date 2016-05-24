update.ep <-
  function(i, y, mn, lis) {
    # update at index \code{i} for the EP approximation given \code{y}
    # (on +1, -1 scale), \code{mn} (on the Gaussian scale) and the current
    # list of parameters \code{lis}.
    
    Sigma <- lis$Sigma
    ttau <- lis$ttau
    tnu <- lis$tnu
    mu <- lis$mu
    
    # calculate approximate cavity parameters \nu_{-i} and \tau_{-i}
    # \sigma^{2}_i = \Sigma_{i, i}
    # \tau_{-i} = \sigma^{-2}_i - \tilde{\tau}_i
    # \nu_{-i} = \mu_i / \sigma^{2}_i + m_i * \tau_{-i} - \tilde{\nu}_i
    sig2_i <- Sigma[i, i]
    tau_ni <- 1 / sig2_i - ttau[i]
    nu_ni <- mu[i] / sig2_i + mn[i] * tau_ni - tnu[i]
    
    # compute marginal moments \hat{\mu}_i and \hat{\sigma}_i^2
    # calculate required derivatives of individual log partition function
    z <- nu_ni / tau_ni / (sqrt(1 + 1 / tau_ni))
    yz <- y[i] * z
    lZ <- pnorm(yz, log.p = TRUE)
    n_p <- dnorm(yz) / exp(lZ)
    dlZ <- y[i] * n_p / sqrt(1 + 1/tau_ni)
    d2lZ <- -n_p * (yz + n_p) / (1 + 1/tau_ni)
    
    # save old \tilde{\tau} before finding new one
    ttau_old <- ttau[i]
    
    # update \tilde{\tau}_i, forcing it non-negative
    # then update \tilde{\nu}_i
    ttau[i] <- max(-d2lZ / (1 + d2lZ / tau_ni), 0)
    tnu[i] <- (dlZ + (mn[i] - nu_ni / tau_ni) * d2lZ) / (1 + d2lZ / tau_ni)
    
    # rank-1 update \Sigma and \mathbf{\mu}
    # \delta\sigma^2
    ds2 <- ttau[i] - ttau_old
    # get column i
    #   si <- Sigma[, i]
    # recompute\Sigma \& \mu
    lis <- update.sigma.mu(Sigma, ds2, i, tnu) 
    
    # return update list of parameters
    return (list(Sigma = lis$Sigma,
                 ttau = ttau,
                 tnu = tnu,
                 mu = lis$mu))
}