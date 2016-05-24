library(unmarked)


sim1 <- function(R=50, J=3, K=3, lambda=5, phi=0.6, p=0.4) {
    M <- rpois(R, lambda) # super-population size
    N <- matrix(NA, R, J) # Population available
    y <- array(NA, c(R, K, J)) # Detected
    for(i in 1:R) {
        for(j in 1:J) {
            N[i,j] <- rbinom(1, M[i], phi)
            y[i,,j] <- rbinom(K, N[i,j], p)
        }
    }
    y <- matrix(y, R)
    return(list(y=y, N=N))
}

set.seed(348)
y1 <- sim1()$y

y1[1,] <- NA
y1[2, 1:3] <- NA
y1[3, 4:6] <- NA
umf <- unmarkedFrameGPC(y=y1, numPrimary=3)

fm1.1 <- gpcount(~1, ~1, ~1, umf, K=40, control=list(trace=TRUE, REPORT=1))
fm1.1r <- gpcount(~1, ~1, ~1, umf, K=40, engine="R",
                  control=list(trace=TRUE, REPORT=1))
fm1.2 <- gpcount(~1, ~1, ~1, umf, K=40, mixture="NB",
                 control=list(trace=TRUE, REPORT=1))


nsim1 <- 10
simout1 <- matrix(NA, nsim1, 3)
lam1 <- 5
phi1 <- 0.5
p1 <- 0.4
nPrimary1 <- 3
set.seed(404)
for(i in 1:nsim1) {
    cat("doing", i, "\n")
    sim1.i <- sim1(lambda=lam1, phi=phi1, p=p1, J=nPrimary1)$y
    umf1.i <- unmarkedFrameGPC(y=sim1.i, numPrimary=nPrimary1)
    fm1.i <- gpcount(~1, ~1, ~1, umf1.i, K=50, engine="C", se=FALSE)
    mle1.i <- coef(fm1.i)
    simout1[i,] <- c(exp(mle1.i[1]), plogis(mle1.i[2:3]))
    cat("  mle =", simout1[i,], "\n")
}

op <- par(mfrow=c(3,1), mai=c(0.5,0.5,0.1,0.1))
hist(simout1[,1]); abline(v=lam1, lwd=2, col=4)
hist(simout1[,2]); abline(v=phi1, lwd=2, col=4)
hist(simout1[,3]); abline(v=p1, lwd=2, col=4)
par(op)













# Covariates

set.seed(568)
R <- 50
J <- 4
K <- 3
x1 <- rnorm(R)
x2 <- matrix(rnorm(R*J), R, J)
x3 <- matrix(rnorm(R*K*J), R, K*J)
x1[2] <- NA
x2[3,] <- NA
x2[4,2] <- NA
x3[5,1:K] <- NA

sim2 <- function(x1, x2, x3,
                 lam0=0, lam1=1, phi0=1, phi1=1, p0=0, p1=1) {
    R <- length(x1)
    J <- ncol(x2)
    K <- ncol(x3)/J
    lambda <- exp(lam0 + lam1*x1)
    phi <- plogis(phi0 + phi1*x2)
    p <- plogis(p0 + p1*x3)
    p <- array(p, c(R, K, J))
    M <- rpois(R, lambda) # super-population size
    N <- matrix(NA, R, J) # Population available
    y <- array(NA, c(R, K, J)) # Detected
    for(i in 1:R) {
        for(j in 1:J) {
            N[i,j] <- rbinom(1, M[i], phi[i,j])
            y[i,,j] <- rbinom(K, N[i,j], p[i,,j])
        }
    }
    y <- matrix(y, R)
    return(list(y=y, N=N))
}

y2 <- sim2(x1, x2, x3)$y
umf2 <- unmarkedFrameGPC(y=y2,
                         siteCovs=data.frame(x1),
                         yearlySiteCovs=list(x2=x2),
                         obsCovs = list(x3=x3), numPrimary=J)
summary(umf2)

fm2.1 <- gpcount(~x1, ~x2, ~x3, umf2, K=40, engine="C",
                 control=list(trace=TRUE, REPORT=1))
fm2.1r <- gpcount(~x1, ~x2, ~x3, umf2, K=40, engine="R",
                 control=list(trace=TRUE, REPORT=1))



nsim2 <- 5
simout2 <- matrix(NA, nsim2, 6)
nPrimary2 <- 4
lam0 <- 0
lam1 <- 1
phi0 <- 1
phi1 <- 1
p0 <- 0
p1 <- 1
set.seed(3434)
for(i in 1:nsim2) {
#    if(i %% 1 == 5)
    cat("doing", i, "\n")
    sim2.i <- sim2(x1, x2, x3, lam0, lam1, phi0, phi1, p0, p1)$y
    umf2.i <- unmarkedFrameGPC(y=sim2.i, siteCovs=data.frame(x1),
                               obsCovs=list(x3=x3),
                               yearlySiteCovs=list(x2=x2),
                               numPrimary=nPrimary2)
    fm2.i <- gpcount(~x1, ~x2, ~x3, umf2.i, K=50, engine="C", se=FALSE)
    mle2.i <- coef(fm2.i)
    simout2[i,] <- mle2.i
    cat("  mle =", mle2.i, "\n")
}

op <- par(mfrow=c(3,2), mai=c(0.5,0.5,0.1,0.1))
hist(simout2[,1]); abline(v=lam0, lwd=2, col=4)
hist(simout2[,2]); abline(v=lam1, lwd=2, col=4)
hist(simout2[,3]); abline(v=phi0, lwd=2, col=4)
hist(simout2[,4]); abline(v=phi1, lwd=2, col=4)
hist(simout2[,5]); abline(v=p0, lwd=2, col=4)
hist(simout2[,6]); abline(v=p1, lwd=2, col=4)
par(op)
