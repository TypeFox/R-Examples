## JAGS model file for the n.opt function

model {
    mu ~ dnorm(nu, 1 / tau^2)
        
    for (i in 1:k) {
        X[i] ~ dnorm(mu, n / sigma^2)
    }    
}
