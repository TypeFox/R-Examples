
library(unmarked)

## Simulate no covariates, constant sampling period intervals,
## no secondary samples

sim1 <- function(lambda=1, gamma=0.5, omega=0.8, p=0.7, M=100, T=5)
{
    y <- N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    for(t in 1:(T-1)) {
        S[,t] <- rbinom(M, N[,t], omega)
        G[,t] <- rpois(M, gamma)
        N[,t+1] <- S[,t] + G[,t]
        }
    y[] <- rbinom(M*T, N, p)
    return(y)
}





set.seed(3223)
nsim1 <- 100
simout1 <- matrix(NA, nsim1, 4)
colnames(simout1) <- c('lambda', 'gamma', 'omega', 'p')
for(i in 1:nsim1) {
    cat("sim1:", i, "\n"); flush.console()
    lambda <- 1
    gamma <- 0.5
    omega <- 0.8
    p <- 0.7
    y.sim1 <- sim1(lambda, gamma, omega, p)
    umf1 <- unmarkedFramePCO(y = y.sim1, numPrimary=5)
    m1 <- pcountOpen(~1, ~1, ~1, ~1, umf1, K=15,
        starts=c(log(lambda), log(gamma), plogis(omega), plogis(p)),
        se=FALSE)
    e <- coef(m1)
    simout1[i, 1:2] <- exp(e[1:2])
    simout1[i, 3:4] <- plogis(e[3:4])
    cat("  mle =", simout1[i,], "\n")
    }

#png("pcountOpenSim1.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout1[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout1[,2], xlab=expression(gamma)); abline(v=gamma, lwd=2, col=4)
hist(simout1[,3], xlab=expression(omega)); abline(v=omega, lwd=2, col=4)
hist(simout1[,4], xlab=expression(p)); abline(v=p, lwd=2, col=4)
#dev.off()





## Simulate covariate model with constant intervals



sim2 <- function(lam=c(0,1), gam=c(-1,-1), om=c(2,-1), p=c(-1,1), M=100,
                 T=5)
{
    y <- gamma <- omega <- det <- N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    veght <- rnorm(M)
    isolation <- matrix(rnorm(M*T), M, T)
    time <- matrix(rnorm(M*T, 1), M, T)
    lambda <- exp(lam[1] + lam[2]*veght)
    gamma[] <- exp(gam[1] + gam[2]*isolation)
    omega[] <- plogis(om[1] + om[2]*isolation)
    det[] <- plogis(p[1] + p[2]*time)

    N[,1] <- rpois(M, lambda)
    for(t in 1:(T-1)) {
        S[,t] <- rbinom(M, N[,t], omega[,t])
        G[,t] <- rpois(M, gamma[,t])
        N[,t+1] <- S[,t] + G[,t]
        }
    y[] <- rbinom(M*T, N, det)
    return(list(y=y, covs=data.frame(veght=veght,
        isolation=isolation, time=time)))
}





nsim2 <- 20
simout2 <- matrix(NA, nsim2, 8)
colnames(simout2) <- c('lam0', 'lam1', 'gam0', 'gam1', 'om0', 'om1',
                       'p0', 'p1')
for(i in 1:nsim2) {
    cat("sim2:", i, "\n")
    lam <- c(-2, 1)
    gam <- c(-1, -1)
    om <- c(0, -1)
    p <- c(-1, 1)
    T <- 5
    sim2out <- sim2(lam, gam, om, p, T=T, M=50)
    y.sim2 <- sim2out$y
    covs <- sim2out$covs
    cn <- colnames(covs)
    siteCovs <- covs[,grep("veght", cn), drop=FALSE]
    yearlySiteCovs <- list(isolation = covs[,grep("isolation", cn)])
    obsCovs <- list(time = covs[,grep("time", cn)])
    umf2 <- unmarkedFramePCO(y = y.sim2, siteCovs=siteCovs,
        yearlySiteCovs=yearlySiteCovs, obsCovs=obsCovs, numPrimary=T)
    m2 <- pcountOpen(~veght, ~isolation, ~isolation, ~time, umf2,
                     K=40, se=FALSE, starts=c(lam, gam, om, p),
                     control=list(trace=TRUE, REPORT=1))
    e <- coef(m2)
    simout2[i, ] <- e
    cat("  mle =", e, "\n")
    }

#png("pcountOpenSim2.png", width=6, height=8, units="in", res=360)
par(mfrow=c(4,2))
hist(simout2[,1], xlab=expression(lambda)); abline(v=lam[1], lwd=2, col=4)
hist(simout2[,2], xlab=expression(lambda)); abline(v=lam[2], lwd=2, col=4)
hist(simout2[,3], xlab=expression(gamma)); abline(v=gam[1], lwd=2, col=4)
hist(simout2[,4], xlab=expression(gamma)); abline(v=gam[2], lwd=2, col=4)
hist(simout2[,5], xlab=expression(omega)); abline(v=om[1], lwd=2, col=4)
hist(simout2[,6], xlab=expression(omega)); abline(v=om[2], lwd=2, col=4)
hist(simout2[,7], xlab=expression(p)); abline(v=p[1], lwd=2, col=4)
hist(simout2[,8], xlab=expression(p)); abline(v=p[2], lwd=2, col=4)
dev.off()













## Simulate uneven sampling period intervals with all dates[i,1]==1

sim3 <- function(lambda=4, gamma=0.1, omega=0.8, p=0.7, M=100, T=5)
{
    y <- N <- date <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    date[,1] <- 1
    for(i in 1:M) {
    for(t in 2:T) {
        delta <- max(rpois(1, 5), 1)
        date[i, t] <- date[i, t-1] + delta
        S[i, t-1] <- rbinom(1, N[i, t-1], omega)
        G[i, t-1] <- rpois(1, gamma)
        N[i, t] <- S[i, t-1] + G[i, t-1]
        if(delta > 1) {
            for(d in 2:delta) {
                S[i, t-1] <- rbinom(1, N[i, t], omega)
                G[i, t-1] <- rpois(1, gamma)
                N[i, t] <- S[i, t-1] + G[i, t-1]
                }
            }
        }}
    y[] <- rbinom(M*T, N, p)
    mode(date) <- "integer"
    return(list(y=y, dates=date))
}









set.seed(373)
nsim3 <- 100
simout3 <- matrix(NA, nsim3, 4)
colnames(simout3) <- c('lambda', 'gamma', 'omega', 'p')
for(i in 1:nsim3) {
    cat("sim3:", i, "\n"); flush.console()
    lambda <- 4
    gamma <- 0.3
    omega <- 0.7
    p <- 0.7
    T <- 5
    yd <- sim3(lambda, gamma, omega, p, M=100, T=5)
    y.sim3 <- yd$y
    dates3 <- yd$dates
    umf3 <- unmarkedFramePCO(y = y.sim3, primaryPeriod=dates3,
                             numPrimary=T)
    m3 <- pcountOpen(~1, ~1, ~1, ~1, umf3, K=20,
        starts=c(log(lambda), log(gamma), plogis(omega), plogis(p)),
                     se=FALSE)
    e <- coef(m3)
    simout3[i, 1:2] <- exp(e[1:2])
    simout3[i, 3:4] <- plogis(e[3:4])
    cat("  mle =", simout3[i,], "\n")
    }



#png("pcountOpenSim3.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout3[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout3[,2], xlab=expression(gamma)); abline(v=gamma, lwd=2, col=4)
hist(simout3[,3], xlab=expression(omega)); abline(v=omega, lwd=2, col=4)
hist(simout3[,4], xlab=expression(p)); abline(v=p, lwd=2, col=4)
dev.off()











# Auto-regressive model




sim4 <- function(lambda=1, gamma=0.5, omega=0.8, p=0.7, M=100, T=5)
{
    y <- N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    for(t in 1:(T-1)) {
        S[,t] <- rbinom(M, N[,t], omega)
        G[,t] <- rpois(M, gamma*N[,t])
        N[,t+1] <- S[,t] + G[,t]
        }
    y[] <- rbinom(M*T, N, p)
    return(y)
}



set.seed(3223)
nsim4 <- 100
simout4 <- matrix(NA, nsim4, 4)
colnames(simout4) <- c('lambda', 'gamma', 'omega', 'p')
for(i in 1:nsim4) {
    cat("sim4:", i, "\n")
    lambda <- 1
    gamma <- 0.5
    omega <- 0.7
    p <- 0.7
    T <- 5
    y.sim4 <- sim4(lambda, gamma, omega, p, T=T)
    umf4 <- unmarkedFramePCO(y = y.sim4, numPrimary=T)
    m4 <- pcountOpen(~1, ~1, ~1, ~1, umf4, K=30, dynamics="autoreg",
                     starts=c(log(lambda), log(gamma), plogis(omega),
                     plogis(p)), se=FALSE)
    e <- coef(m4)
    simout4[i, 1:2] <- exp(e[1:2])
    simout4[i, 3:4] <- plogis(e[3:4])
    cat("  mle =", simout4[i,], "\n")
    }

#png("pcountOpenSim4.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout4[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout4[,2], xlab=expression(gamma)); abline(v=gamma, lwd=2, col=4)
hist(simout4[,3], xlab=expression(omega)); abline(v=omega, lwd=2, col=4)
hist(simout4[,4], xlab=expression(p)); abline(v=p, lwd=2, col=4)
dev.off()










# No trend model




sim5 <- function(lambda=1, omega=0.8, p=0.7, M=100, T=5)
{
    y <- N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    gamma <- (1-omega)*lambda
    for(t in 1:(T-1)) {
        S[,t] <- rbinom(M, N[,t], omega)
        G[,t] <- rpois(M, gamma)
        N[,t+1] <- S[,t] + G[,t]
        }
    y[] <- rbinom(M*T, N, p)
    return(y)
}



set.seed(3223)
nsim5 <- 100
simout5 <- matrix(NA, nsim5, 3)
colnames(simout5) <- c('lambda', 'omega', 'p')
for(i in 1:nsim5) {
    cat("sim5:", i, "\n"); flush.console()
    lambda <- 1
    omega <- 0.7
    p <- 0.7
    T <- 5
    y.sim5 <- sim5(lambda, omega, p, T=T)
    umf5 <- unmarkedFramePCO(y = y.sim5, numPrimary=T)
    m5 <- pcountOpen(~1, ~1, ~1, ~1, umf5, K=20, dynamics="notrend",
        starts=c(log(lambda), plogis(omega), plogis(p)), se=FALSE)
    e <- coef(m5)
    simout5[i, 1] <- exp(e[1])
    simout5[i, 2:3] <- plogis(e[2:3])
    cat("  mle =", simout5[i,], "\n")
    }

#png("pcountOpenSim5.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout5[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout5[,2], xlab=expression(omega)); abline(v=omega, lwd=2, col=4)
hist(simout5[,3], xlab=expression(p)); abline(v=p, lwd=2, col=4)
#dev.off()












## Simulate data with some dates[i,1] > 1


sim6 <- function(lambda=4, gamma=0.1, omega=0.8, p=0.7, M=100, T=5)
{
    y <- N <- date <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    date[,1] <- pmax(rpois(M, 5), 1)

    for(i in 1:M) {
    if(date[i,1] > 1) {
        for(d in 2:date[i, 1]) {
            S1 <- rbinom(1, N[i, 1], omega)
            G1 <- rpois(1, gamma)
            N[i, 1] <- S1 + G1
            }
        }
    for(t in 2:T) {
        date[i, t] <- date[i, t-1] + 1
        S[i, t-1] <- rbinom(1, N[i, t-1], omega)
        G[i, t-1] <- rpois(1, gamma)
        N[i, t] <- S[i, t-1] + G[i, t-1]
        }}
    y[] <- rbinom(M*T, N, p)
    mode(date) <- "integer"
    return(list(y=y, dates=date))
}




set.seed(3223)
nsim6 <- 100
simout6 <- matrix(NA, nsim6, 4)
colnames(simout6) <- c('lambda', 'gamma', 'omega', 'p')
for(i in 1:nsim6) {
    cat("sim6:", i, "\n")
    lambda <- 1
    gamma <- 0.5
    omega <- 0.8
    p <- 0.7
    T <- 5
    yd <- sim6(lambda, gamma, omega, p, M=100, T=T)
    y.sim6 <- yd$y
    dates6 <- yd$dates
    umf6 <- unmarkedFramePCO(y = y.sim6, primaryPeriod=dates6,
                             numPrimary=T)
    m6 <- pcountOpen(~1, ~1, ~1, ~1, umf6, K=25,
        starts=c(log(lambda), log(gamma), plogis(omega), plogis(p)),
                     se=FALSE)
    e <- coef(m6)
    simout6[i, 1:2] <- exp(e[1:2])
    simout6[i, 3:4] <- plogis(e[3:4])
    cat("  mle =", simout6[i,], "\n")
    }



#png("pcountOpenSim6.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout6[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout6[,2], xlab=expression(gamma)); abline(v=gamma, lwd=2, col=4)
hist(simout6[,3], xlab=expression(omega)); abline(v=omega, lwd=2, col=4)
hist(simout6[,4], xlab=expression(p)); abline(v=p, lwd=2, col=4)
#dev.off()














## Simulate no covariates, constant sampling period intervals,
## WITH secondary samples

sim7 <- function(lambda=1, gamma=0.5, omega=0.8, p=0.7, M=100, T=5, J=3)
{
    y <- matrix(NA, M, J*T)
    N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    for(t in 1:(T-1)) {
        S[,t] <- rbinom(M, N[,t], omega)
        G[,t] <- rpois(M, gamma)
        N[,t+1] <- S[,t] + G[,t]
        }
    N <- N[,rep(1:T, each=J)]
    y[] <- rbinom(M*J*T, N, p)
    return(y)
}



library(unmarked)
set.seed(3223)
nsim7 <- 100
simout7 <- matrix(NA, nsim7, 4)
colnames(simout7) <- c('lambda', 'gamma', 'omega', 'p')
for(i in 1:nsim7) {
    cat("sim7:", i, "\n")
    lambda <- 1
    gamma <- 0.5
    omega <- 0.8
    p <- 0.7
    T <- 5
    y.sim7 <- sim7(lambda, gamma, omega, p, T=T)
    umf7 <- unmarkedFramePCO(y = y.sim7, numPrimary=T)
    m7 <- pcountOpen(~1, ~1, ~1, ~1, umf7, K=15,
              starts=c(log(lambda), log(gamma), plogis(omega), plogis(p)),
              se=FALSE)
    e <- coef(m7)
    simout7[i, 1:2] <- exp(e[1:2])
    simout7[i, 3:4] <- plogis(e[3:4])
    cat("mle = ", simout7[i,], "\n")
    }

#png("pcountOpenSim7.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout7[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout7[,2], xlab=expression(gamma)); abline(v=gamma, lwd=2, col=4)
hist(simout7[,3], xlab=expression(omega)); abline(v=omega, lwd=2, col=4)
hist(simout7[,4], xlab=expression(p)); abline(v=p, lwd=2, col=4)
#dev.off()















## Simulate covariate model with constant intervals and secondary samples



sim8 <- function(lam=c(0,1), gam=c(-1,-1), om=c(2,-1), p=c(-1,1), M=100,
    T=5, J=3)
{
    y <- det <- matrix(NA, M, J*T)
    gamma <- omega <- N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    veght <- rnorm(M)
    isolation <- matrix(rnorm(M*T), M, T)
    time <- matrix(rnorm(M*J*T, 1), M, J*T)
    lambda <- exp(lam[1] + lam[2]*veght)
    gamma[] <- exp(gam[1] + gam[2]*isolation)
    omega[] <- plogis(om[1] + om[2]*isolation)
    det[] <- plogis(p[1] + p[2]*time)

    N[,1] <- rpois(M, lambda)
    for(t in 1:(T-1)) {
        S[,t] <- rbinom(M, N[,t], omega[,t])
        G[,t] <- rpois(M, gamma[,t])
        N[,t+1] <- S[,t] + G[,t]
        }
    N <- N[,rep(1:T, each=J)]
    y[] <- rbinom(M*J*T, N, det)
    return(list(y=y, covs=data.frame(veght=veght,
        isolation=isolation, time=time)))
}





nsim8 <- 100
simout8 <- matrix(NA, nsim8, 8)
colnames(simout8) <- c('lam0', 'lam1', 'gam0', 'gam1', 'om0', 'om1', 'p0', 'p1')
for(i in 1:nsim8) {
    cat("sim8:", i, "\n"); flush.console()
    lam <- c(-2, 1)
    gam <- c(-1, -1)
    om <- c(0, -1)
    p <- c(-1, 1)
    T <- 5
    sim8out <- sim8(lam, gam, om, p, T=T)
    y.sim8 <- sim8out$y
    covs <- sim8out$covs
    cn <- colnames(covs)
    siteCovs <- covs[,grep("veght", cn), drop=FALSE]
    yearlySiteCovs <- list(isolation=covs[,grep("isolation", cn)])
    obsCovs <- list(time = covs[,grep("time", cn)])
    umf8 <- unmarkedFramePCO(y = y.sim8, siteCovs=siteCovs,
        yearlySiteCovs=yearlySiteCovs, obsCovs=obsCovs, numPrimary=T)
    m8 <- pcountOpen(~veght, ~isolation, ~isolation, ~time, umf8, K=30,
                     se=F, starts=c(lam, gam, om, p))
    e <- coef(m8)
    simout8[i, ] <- e
    cat("  mle=", e, "\n")
    }

#png("pcountOpenSim8.png", width=6, height=8, units="in", res=360)
par(mfrow=c(4,2))
hist(simout8[,1], xlab=expression(lambda)); abline(v=lam[1], lwd=2, col=4)
hist(simout8[,2], xlab=expression(lambda)); abline(v=lam[2], lwd=2, col=4)
hist(simout8[,3], xlab=expression(gamma)); abline(v=gam[1], lwd=2, col=4)
hist(simout8[,4], xlab=expression(gamma)); abline(v=gam[2], lwd=2, col=4)
hist(simout8[,5], xlab=expression(omega)); abline(v=om[1], lwd=2, col=4)
hist(simout8[,6], xlab=expression(omega)); abline(v=om[2], lwd=2, col=4)
hist(simout8[,7], xlab=expression(p)); abline(v=p[1], lwd=2, col=4)
hist(simout8[,8], xlab=expression(p)); abline(v=p[2], lwd=2, col=4)
#dev.off()















## Simulate no covariates, constant sampling period intervals,
## WITH secondary samples and missing values

sim9 <- function(lambda=1, gamma=0.5, omega=0.8, p=0.7, M=100, T=5, J=3,
                 nMissing=50)
{
    y <- matrix(NA, M, J*T)
    N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    for(t in 1:(T-1)) {
        S[,t] <- rbinom(M, N[,t], omega)
        G[,t] <- rpois(M, gamma)
        N[,t+1] <- S[,t] + G[,t]
        }
    N <- N[,rep(1:T, each=J)]
    y[] <- rbinom(M*J*T, N, p)
    y[sample.int(M*J*T, nMissing)] <- NA
    return(y)
}



library(unmarked)
set.seed(3223)
nsim9 <- 100
simout9 <- matrix(NA, nsim9, 4)
colnames(simout9) <- c('lambda', 'gamma', 'omega', 'p')
for(i in 1:nsim9) {
    cat("sim9:", i, "\n")
    lambda <- 1
    gamma <- 0.5
    omega <- 0.8
    p <- 0.7
    T <- 5
    y.sim9 <- sim9(lambda, gamma, omega, p, T=T, nMissing=100)
    umf9 <- unmarkedFramePCO(y = y.sim9, numPrimary=T)
    m9 <- pcountOpen(~1, ~1, ~1, ~1, umf9, K=15,
              starts=c(log(lambda), log(gamma), plogis(omega), plogis(p)),
              se=FALSE)
    e <- coef(m9)
    simout9[i, 1:2] <- exp(e[1:2])
    simout9[i, 3:4] <- plogis(e[3:4])
    cat("mle = ", simout9[i,], "\n")
    }

#png("pcountOpenSim9.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout9[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout9[,2], xlab=expression(gamma)); abline(v=gamma, lwd=2, col=4)
hist(simout9[,3], xlab=expression(omega)); abline(v=omega, lwd=2, col=4)
hist(simout9[,4], xlab=expression(p)); abline(v=p, lwd=2, col=4)
#dev.off()

















## Simulate covariate model with secondary samples and missing values



sim10 <- function(lam=c(0,1), gam=c(-1,-1), om=c(2,-1), p=c(-1,1), M=100,
    T=5, J=3, nMissing=50)
{
    y <- det <- matrix(NA, M, J*T)
    gamma <- omega <- N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    veght <- rnorm(M)
    isolation <- matrix(rnorm(M*T), M, T)
    time <- matrix(rnorm(M*J*T, 1), M, J*T)
    lambda <- exp(lam[1] + lam[2]*veght)
    gamma[] <- exp(gam[1] + gam[2]*isolation)
    omega[] <- plogis(om[1] + om[2]*isolation)
    det[] <- plogis(p[1] + p[2]*time)

    N[,1] <- rpois(M, lambda)
    for(t in 1:(T-1)) {
        S[,t] <- rbinom(M, N[,t], omega[,t])
        G[,t] <- rpois(M, gamma[,t])
        N[,t+1] <- S[,t] + G[,t]
        }
    N <- N[,rep(1:T, each=J)]
    y[] <- rbinom(M*J*T, N, det)
    na.ind <- sample.int(M*J*T, nMissing)
    veght[na.ind<=M] <- NA
    isolation[na.ind<=(M*T)] <- NA
    time[na.ind] <- NA
    covs <- data.frame(veght=veght, isolation=isolation, time=time)
    y[na.ind] <- NA
    return(list(y=y, covs=covs))
}




set.seed(4483499)
nsim10 <- 100
simout10 <- matrix(NA, nsim10, 7)
colnames(simout10) <- c('lam0', 'lam1', 'gam0', 'gam1', 'om0', #'om1',
                        'p0', 'p1')
for(i in 1:nsim10) {
    cat("sim10", i, "\n"); flush.console()
    lam <- c(-2, 1)
    gam <- c(-1, -1)
    om <- c(0, 0)
    p <- c(-1, 1)
    T <- 5
    sim10out <- sim10(lam, gam, om, p, T=T)
    y.sim10 <- sim10out$y
    covs <- sim10out$covs
    cn <- colnames(covs)
    siteCovs <- covs[,grep("veght", cn), drop=FALSE]
    yearlySiteCovs <- list(isolation=covs[,grep("isolation", cn)])
    obsCovs <- list(time = covs[,grep("time", cn)])
    umf10 <- unmarkedFramePCO(y = y.sim10, siteCovs=siteCovs,
        yearlySiteCovs=yearlySiteCovs, obsCovs=obsCovs, numPrimary=T)
    m10 <- pcountOpen(~veght, ~isolation, ~1, ~time, umf10, K=30,
                      se=F, starts=c(lam, gam, 0, p),
                      control=list(trace=F, REPORT=1))
    e <- coef(m10)
    simout10[i, ] <- e
    cat("  mle=", e, "\n")
    }

#png("pcountOpenSim10.png", width=6, height=10, units="in", res=360)
par(mfrow=c(4,2))
hist(simout10[,1], xlab=expression(lambda)); abline(v=lam[1], lwd=2, col=4)
hist(simout10[,2], xlab=expression(lambda)); abline(v=lam[2], lwd=2, col=4)
hist(simout10[,3], xlab=expression(gamma)); abline(v=gam[1], lwd=2, col=4)
hist(simout10[,4], xlab=expression(gamma)); abline(v=gam[2], lwd=2, col=4)
hist(simout10[,5], xlab=expression(omega)); abline(v=om[1], lwd=2, col=4)
hist(simout10[,6], xlab=expression(omega)); abline(v=om[2], lwd=2, col=4)
hist(simout10[,7], xlab=expression(p)); abline(v=p[1], lwd=2, col=4)
hist(simout10[,8], xlab=expression(p)); abline(v=p[2], lwd=2, col=4)
#dev.off()


m10 <- pcountOpen(~veght, ~1, ~1, ~time, umf10, K=30,
                  se=F, #starts=c(0, 0, 0, p),
                  control=list(trace=TRUE, REPORT=1))




trace(unmarked:::handleNA, browser, browser, signature="unmarkedFramePCO")
untrace(unmarked:::handleNA, signature="unmarkedFramePCO")


debugonce(pcountOpen)



## Simulate "trend model", no covariates, constant intervals,
## no secondary samples

sim11 <- function(lambda=1, gamma=0.5, p=0.7, M=100, T=5)
{
    y <- N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    for(t in 2:T) {
        N[,t] <- rpois(M, gamma*N[,t-1])
        }
    y[] <- rbinom(M*T, N, p)
    return(y)
}





set.seed(3223)
nsim11 <- 100
simout11 <- matrix(NA, nsim11, 3)
colnames(simout11) <- c('lambda', 'gamma', 'p')
for(i in 1:nsim11) {
    cat("sim11:", i, "\n"); flush.console()
    lambda <- 2
    gamma <- 0.5
    p <- 0.7
    y.sim11 <- sim11(lambda, gamma, p)
    umf11 <- unmarkedFramePCO(y = y.sim11, numPrimary=5)
    m11 <- pcountOpen(~1, ~1, ~1, ~1, umf11, K=40, dynamics="trend",
        starts=c(log(lambda), log(gamma), plogis(p)),
        se=FALSE)
    e <- coef(m11)
    simout11[i, 1:2] <- exp(e[1:2])
    simout11[i, 3] <- plogis(e[3])
    cat("  mle =", simout11[i,], "\n")
    }

#png("pcountOpenSim1.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout11[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout11[,2], xlab=expression(gamma)); abline(v=gamma, lwd=2, col=4)
hist(simout11[,3], xlab=expression(p)); abline(v=p, lwd=2, col=4)
#dev.off()











## Simulate trend model with ZIP dist

sim12 <- function(lambda=1, gamma=0.5, p=0.7, psi=0.3, M=100, T=5)
{
    y <- N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    N[runif(M) < psi, 1] <- 0
    for(t in 2:T) {
        N[,t] <- rpois(M, gamma*N[,t-1])
        }
    y[] <- rbinom(M*T, N, p)
    return(y)
}





set.seed(3223)
nsim12 <- 100
simout12 <- matrix(NA, nsim12, 4)
colnames(simout12) <- c('lambda', 'gamma', 'p', 'psi')
for(i in 1:nsim12) {
    cat("sim12:", i, "\n")
    lambda <- 2
    gamma <- 0.5
    p <- 0.7
    psi <- 0.3
    y.sim12 <- sim12(lambda, gamma, p, psi)
    umf12 <- unmarkedFramePCO(y = y.sim12, numPrimary=5)
    m12 <- pcountOpen(~1, ~1, ~1, ~1, umf12, K=40, dynamics="trend",
                      mixture="ZIP",
        starts=c(log(lambda), log(gamma), plogis(p), plogis(psi)),
        se=FALSE)
    e <- coef(m12)
    simout12[i, 1:2] <- exp(e[1:2])
    simout12[i, 3:4] <- plogis(e[3:4])
    cat("  mle =", simout12[i,], "\n")
    }

#png("pcountOpenSim1.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout12[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout12[,2], xlab=expression(gamma)); abline(v=gamma, lwd=2, col=4)
hist(simout12[,3], xlab=expression(p)); abline(v=p, lwd=2, col=4)
hist(simout12[,4], xlab=expression(psi)); abline(v=psi, lwd=2, col=4)
#dev.off()



## Simulate Ricker model 

sim13 <- function(lambda=1, gamma=0.1, omega=1.5, p=0.7, M=100, T=5)
{
    y <- N <- matrix(NA, M, T)
    N[,1] <- rpois(M, lambda)
    for(t in 2:T) {
        N[,t] <- rpois(M, N[,t-1]*exp(gamma*(1-N[,t-1]/omega)))
    }
    y[] <- rbinom(M*T, N, p)
    return(y)
}





set.seed(3223)
nsim13 <- 100
simout13 <- matrix(NA, nsim13, 4)
colnames(simout13) <- c('lambda', 'gamma', 'omega', 'p')
for(i in 1:nsim13) {
    cat("sim13:", i, "\n")
    lambda <- 2
    gamma <- 0.25
    omega <- 2.3
    p <- 0.7
    y.sim13 <- sim13(lambda, gamma, omega, p)
    umf13 <- unmarkedFramePCO(y = y.sim13, numPrimary=5)
    m13 <- pcountOpen(~1, ~1, ~1, ~1, umf13, K=40, dynamics="ricker",
        starts=c(log(lambda), log(gamma), log(omega), plogis(p)),
        se=FALSE)
    e <- coef(m13)
    simout13[i, 1:3] <- exp(e[1:3])
    simout13[i, 4] <- plogis(e[4])
    cat("  mle =", simout13[i,], "\n")
    }

#png("pcountOpenSim1.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout13[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout13[,2], xlab=expression(gamma)); abline(v=gamma, lwd=2, col=4)
hist(simout13[,3], xlab=expression(omega)); abline(v=omega, lwd=2, col=4)
hist(simout13[,4], xlab=expression(p)); abline(v=p, lwd=2, col=4)
#dev.off()



## Simulate Gompertz model 

sim14 <- function(lambda=1, gamma=0.1, omega=1.5, p=0.7, M=100, T=5)
{
    y <- N <- matrix(NA, M, T)
    N[,1] <- rpois(M, lambda)
    for(t in 2:T) {
        N[,t] <- rpois(M, N[,t-1]*exp(gamma*(1-log(N[,t-1]+1)/log(omega+1))))
    }
    y[] <- rbinom(M*T, N, p)
    return(y)
}

set.seed(3223)
nsim14 <- 100
simout14 <- matrix(NA, nsim14, 4)
colnames(simout14) <- c('lambda', 'gamma', 'omega', 'p')
for(i in 1:nsim14) {
    cat("sim14:", i, "\n")
    lambda <- 2
    gamma <- 0.25
    omega <- 2.3
    p <- 0.7
    y.sim14 <- sim14(lambda, gamma, omega, p)
    umf14 <- unmarkedFramePCO(y = y.sim14, numPrimary=5)
    m14 <- pcountOpen(~1, ~1, ~1, ~1, umf14, K=40, dynamics="gompertz",
        starts=c(log(lambda), log(gamma), log(omega), plogis(p)),
        se=FALSE)
    e <- coef(m14)
    simout14[i, 1:3] <- exp(e[1:3])
    simout14[i, 4] <- plogis(e[4])
    cat("  mle =", simout14[i,], "\n")
    }

#png("pcountOpenSim1.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout14[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout14[,2], xlab=expression(gamma)); abline(v=gamma, lwd=2, col=4)
hist(simout14[,3], xlab=expression(omega)); abline(v=omega, lwd=2, col=4)
hist(simout14[,4], xlab=expression(p)); abline(v=p, lwd=2, col=4)
#dev.off()



## Simulate trend + immigration model 

sim15 <- function(lambda=1, gamma=0.5, iota=1, p=0.7, M=100, T=5)
{
    y <- N <- matrix(NA, M, T)
    N[,1] <- rpois(M, lambda)
    for(t in 2:T) {
        N[,t] <- rpois(M, gamma*N[,t-1] + iota)
        }
    y[] <- rbinom(M*T, N, p)
    return(y)
}

set.seed(3223)
nsim15 <- 100
simout15 <- matrix(NA, nsim15, 4)
colnames(simout15) <- c('lambda', 'gamma', 'iota', 'p')
for(i in 1:nsim15) {
    cat("sim15:", i, "\n")
    lambda <- 2
    gamma <- 0.25
    iota <- 0.5
    p <- 0.7
    y.sim15 <- sim15(lambda, gamma, iota, p)
    umf15 <- unmarkedFramePCO(y = y.sim15, numPrimary=5)
    m15 <- pcountOpen(~1, ~1, ~1, ~1, umf15, K=40, dynamics="trend",
        starts=c(log(lambda), log(gamma), plogis(p), log(iota)),
        se=TRUE, immigration=TRUE, iotaformula=~1)
    e <- coef(m15)
    simout15[i, 1:3] <- exp(e[c(1:2,4)])
    simout15[i, 4] <- plogis(e[3])
    cat("  mle =", simout15[i,], "\n")
    }

#png("pcountOpenSim1.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout15[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout15[,2], xlab=expression(gamma)); abline(v=gamma, lwd=2, col=4)
hist(simout15[,3], xlab=expression(iota)); abline(v=iota, lwd=2, col=4)
hist(simout15[,4], xlab=expression(p)); abline(v=p, lwd=2, col=4)
dev.off()








