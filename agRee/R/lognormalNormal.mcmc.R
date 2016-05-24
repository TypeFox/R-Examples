lognormalNormal.mcmc <- function(X, nmcmc){

    nsubj <- nrow(X)
    nrater <- ncol(X)
    dd <- data.frame(reading=as.vector(X),
                     method=factor(rep(1:nrater, each=nsubj)),
                     subject=factor(rep(1:nsubj, nrater)))
    form.lmer <- lFormula(reading ~ method + (1|subject), data=dd)
    fit.lmer <- lmer(reading ~ method + (1|subject), data=dd)
                   

    # Set up data
    X <- form.lmer$X
    X <- X[1:nrow(X), 1:ncol(X)]
    P <- ncol(X)

    Z  <- t(as.matrix(form.lmer$reTrms$Zt))
    K <- ncol(Z)

    reading <- dd$reading
    N <- length(reading)

    data=list("reading","X", "Z", "N","K","P")

    # Obtain initial values
    dd <- matrix(dd$reading, ncol=nsubj, byrow=TRUE)
    sRater <- rowMeans(dd)
    betajs <- sRater[2:nrater] - sRater[1]
    vv <- VarCorr(fit.lmer)
    sigma2.e  <- (attributes(vv)$sc)^2

    est <- dd[1,] - residuals(fit.lmer)[1:nsubj]
    beta0 <- min(est)
    alpha.est <- est - beta0
    mu.alpha <- mean(log(alpha.est + 1))
    sigma2.alpha <- var(log(alpha.est + 1))

    inits = function(){list(beta=c(beta0, betajs),
                            b=alpha.est+1,
                            taueps=1/sigma2.e,
                            mub=mu.alpha, taub=1/sigma2.alpha)}
     
    # Fit the model
    cpath <- getwd()
    setwd(cpath)
    cat(
    "model {
        for(i in 1:N){
            reading[i] ~ dnorm(m[i], taueps)
            m[i] <- inprod(beta[], X[i,])+inprod(b[], Z[i,])
        }
        for(k in 1:K){
            b[k] ~ dlnorm(mub, taub)
        }
        for(j in 1:P){
            beta[j] ~ dnorm(0, 1.0E-6)
        }
        taueps ~ dgamma(1.0E-3, 1.0E-3)
        mub ~ dnorm(0, 1.0E-6)
        taub ~ dgamma(1.0E-3, 1.0E-3) 
     }",
        file="lognormal-normal-lme.jags")
    
    jagsfit <- jags(data, inits,
                    parameters.to.save=c("beta", "taueps", "mub", "taub"), 
                    model.file="lognormal-normal-lme.jags",
                    n.chains = 1, n.iter=nmcmc, n.burnin=nmcmc/10, n.thin=10)
    unlink("lognormal-normal-lme.jags")
    
    # Extract parameters and obtain the estimates
    mu <- jagsfit$BUGSoutput$sims.list$mub
    sigma2.alpha <- 1 / jagsfit$BUGSoutput$sims.list$taub
    sigma2.e <- 1 / jagsfit$BUGSoutput$sims.list$taueps
    betas <- cbind(jagsfit$BUGSoutput$sims.list$beta[,1],
                   jagsfit$BUGSoutput$sims.list$beta[,2:nrater]+
                   jagsfit$BUGSoutput$sims.list$beta[,1])

    list(mu=mu, sigma2.alpha=sigma2.alpha, sigma2.e=sigma2.e, betas=betas)
    
}
