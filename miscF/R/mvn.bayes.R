 
simMvnConjugate <- function(X, nsim){
    
    n <- nrow(X)
    k <- ncol(X)
    
    est <- mvn.ub(X)
    k0 <- 1
    bar.x <- est$hatMu
    mu0 <- bar.x
    v0 <- k
    Gamma0 <- est$hatSigma
    S <- est$hatSigma*(n-1)

    kn <- k0 + n
    mun <- (k0/kn) * mu0  +  (n/kn) * bar.x
    vn <- v0 + n
    Gamman <- Gamma0 + S + (k0*n / kn) * (bar.x-mu0) %*% t(bar.x-mu0)

    Sigma <- array(0, dim=c(k, k, nsim))
    Mu <- matrix(0, ncol=k, nrow=nsim)

    
    for(i in 1:nsim){
        Sigma[,,i] <- riwish(vn, Gamman)
        Mu[i,] <- mvrnorm(1, mun, Sigma[,,i]/kn)
    }
    list(Sigma=Sigma, Mu=Mu)
}
#----------------------------------------------------------------------

simBvnJeffreys <- function(X, nsim){
    a <- 1
    b <- 0

    k <- ncol(X)
    n <- nrow(X)
    bar.x <- colMeans(X)
    bar.x1 <- bar.x[1] 
    bar.x2 <- bar.x[2]
    S <- rowSums(apply(X, 1, function(x) (x - bar.x) %*% t(x - bar.x)))
    S <- matrix(S, nrow=2, ncol=2)
    s11 <- S[1,1]
    s22 <- S[2,2]
    r <- S[1,2]/sqrt(s11*s22)

   
    Z1 <- rnorm(nsim)
    Z2 <- rnorm(nsim)
    Z3 <- rnorm(nsim)
    chisq.na <- rchisq(nsim, n-a)
    chisq.nb <- rchisq(nsim, n-b)
 
    mu1 <- bar.x1 + 
           Z1/sqrt(chisq.na) * sqrt(s11/n)

    mu2 <- bar.x2 + 
           Z1/sqrt(chisq.na) * r*sqrt(s22/n) + 
           ( Z2/sqrt(chisq.nb) - 
             Z3/sqrt(chisq.nb) * Z1/sqrt(chisq.na)) * 
             sqrt(s22*(1-r^2)/n)

    sigma1 <- sqrt(s11/chisq.na)

    sigma2 <- sqrt(s22*(1-r^2)) *
              sqrt(1/chisq.nb + 
                   1/chisq.na * (Z3/sqrt(chisq.nb) - r/sqrt(1-r^2))^2)
   
    Y <- -Z3/sqrt(chisq.na) + 
             sqrt(chisq.nb)/sqrt(chisq.na) * r/sqrt(1-r^2)
    rho <- Y / sqrt(1+Y^2)

    Mu <- cbind(mu1, mu2)
    Sigma <- array(0, dim=c(k,k,nsim))
    for(i in 1:nsim){
        Sigma[,,i] <- matrix(c(sigma1[i]^2, sigma1[i]*sigma2[i]*rho[i],
                               sigma1[i]*sigma2[i]*rho[i], sigma2[i]^2), 2, 2)
    }
    
    list(Sigma=Sigma, Mu=Mu)
}

simMvnJeffreys <- function(X, nsim){
    n <- nrow(X)
    k <- ncol(X)
    bar.x <- colMeans(X)
    S <- rowSums(apply(X, 1, function(x) (x - bar.x) %*% t(x - bar.x)))
    S <- matrix(S, nrow=k, ncol=k)

    Sigma <- array(0, dim=c(k, k, nsim))
    Mu <- matrix(0, ncol=k, nrow=nsim)
    for(i in 1:nsim){
        Sigma[,,i] <- riwish(n, S)
        Mu[i,] <- mvrnorm(1, bar.x, Sigma[,,i]/n)
    }
    list(Sigma=Sigma, Mu=Mu)
}

mvn.bayes <- function(X, nsim, prior=c("Jeffreys", "Conjugate")){

    k <- ncol(X)
    prior <- match.arg(prior)
    if(prior=="Jeffreys"){
        if(k==2){
            mcmc.results <- simBvnJeffreys(X, nsim)
        }
        else{
            mcmc.results <- simMvnJeffreys(X, nsim)
        }
    }
    else{
        mcmc.results <- simMvnConjugate(X, nsim)
    }
    list(Mu.save=mcmc.results$Mu, Sigma.save=mcmc.results$Sigma)
}
