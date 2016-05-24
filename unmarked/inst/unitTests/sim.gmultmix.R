
# -------------------------- Null Poisson removal model ------------------

set.seed(26)

n <- 50  # number of sites
T <- 4    # number of primary periods
J <- 3    # number of secondary periods

lam <- 3
phi <- 0.5
p <- 0.3

y <- array(NA, c(n, T, J))
M <- rpois(n, lam)          # Local population size
N <- matrix(NA, n, T)       # Individuals availabe for detection

for(i in 1:n) {
    N[i,] <- rbinom(T, M[i], phi)
    y[i,,1] <- rbinom(T, N[i,], p)
    Nleft1 <- N[i,] - y[i,,1]
    y[i,,2] <- rbinom(T, Nleft1, p)
    Nleft2 <- Nleft1 - y[i,,2]
    y[i,,3] <- rbinom(T, Nleft2, p)
    }

y.ijt <- cbind(y[,1,], y[,2,], y[,3,], y[,4,])
umf1 <- unmarkedFrameGMM(y=y.ijt, numPrimary=T, type="removal")


system.time(m1 <- gmultmix(~1, ~1, ~1, data=umf1)) #2.3

# Test 1
checkEqualsNumeric(coef(m1), c(1.3923561, -0.3183231, -0.7864098),
    tolerance=1e-5)

SSE(m1)

(pb1 <- parboot(m1, nsim=50, report=5))
plot(pb1)




# -------------------------- Null NegBin removal model -------------------

set.seed(73)

n <- 50  # number of sites
T <- 4    # number of primary periods
J <- 3    # number of secondary periods

lam <- 3
phi <- 0.5
p <- 0.3
alpha <- 2

y <- array(NA, c(n, T, J))
M <- rnbinom(n, mu=lam, size=alpha)   # Local population size
N <- matrix(NA, n, T)               # Individuals availabe for detection

for(i in 1:n) {
    N[i,] <- rbinom(T, M[i], phi)
    y[i,,1] <- rbinom(T, N[i,], p)
    Nleft1 <- N[i,] - y[i,,1]
    y[i,,2] <- rbinom(T, Nleft1, p)
    Nleft2 <- Nleft1 - y[i,,2]
    y[i,,3] <- rbinom(T, Nleft2, p)
    }

y.ijt <- cbind(y[,1,], y[,2,], y[,3,], y[,4,])
umf2 <- unmarkedFrameGMM(y=y.ijt, numPrimary=T, type="removal")

system.time(m2 <- gmultmix(~1, ~1, ~1, data=umf2, mixture="NB")) #2.3

backTransform(m2, type="alpha")

# Test
checkEqualsNumeric(coef(m2), c(1.118504, 1.414340, -1.394736, 1.056084),
    tol=1e-5)

(pb2 <- parboot(m2, nsim=50, report=5))
plot(pb2)








# --------------------- Poisson removal model w/ covariates --------------

set.seed(37)

n <- 50   # number of sites
T <- 4    # number of primary periods
J <- 3    # number of secondary periods

sc <- rnorm(n)
ysc <- rnorm(n*T)
ysc <- matrix(ysc, n, T)
yr <- factor(rep(1:T, n))
oc <- rnorm(n*J*T)
oc <- array(oc, c(n, J, T))
int <- matrix(1:(T*J), nrow=n, ncol=T*J, byrow=TRUE)
pi <- array(NA, c(n, J, T))

lam <- exp(-1 + 1*sc)
phi <- plogis(2 + -2*ysc)
p <- plogis(1 + -1*oc)

y <- array(NA, c(n,J,T))
M <- rpois(n, lam)
N <- matrix(NA, n, T)

for(i in 1:n) {
    N[i,] <- rbinom(T, M[i], phi[i,])
    y[i,1,] <- rbinom(T, N[i,], p[i,1,])
    Nleft1 <- N[i,] - y[i,1,]
    y[i,2,] <- rbinom(T, Nleft1, p[i,2,])
    Nleft2 <- Nleft1 - y[i,2,]
    y[i,3,] <- rbinom(T, Nleft2, p[i,3,])
    }

umf3 <- unmarkedFrameGMM(y=matrix(y, nrow=n),
    siteCovs = data.frame(sc=sc),
    obsCovs=list(oc=matrix(oc, nrow=n), int=int),
    yearlySiteCovs=data.frame(ysc=as.numeric(t(ysc)), yr=yr),
    numPrimary=T, type="removal")

(m3 <- gmultmix(~sc, ~ysc, ~oc, umf3))
#system.time(m3 <- gmultmix(~sc, ~ysc, ~oc, umf3)) # 4.8

# Test
checkEqualsNumeric(coef(m3), c(-1.2513974, 1.3585940, 2.2889517, -2.1197854,
    1.0450782, -0.8627125), tol=1e-5)

(pb3 <- parboot(m3, nsim=50, report=5))




umf4 <- unmarkedFrameGMM(y=matrix(y, nrow=n),
    siteCovs = data.frame(sc=sc),
    obsCovs=list(oc=matrix(oc, nrow=n), int=int),
    yearlySiteCovs=list(ysc=ysc),
    numPrimary=T, type="removal")










# ------------------------- independent double observer ------------------



sim.doub <- function(nSites=200, nReps=2, lambda=1, phi=0.6,
                     pA=0.8, pB=0.6, alpha=0.5)
{

    N <- matrix(NA, nSites, nReps)
    y <- array(NA, c(nSites, 3, nReps))

    # Abundance at each site (quadrat)
    M <- rnbinom(nSites, size=alpha, mu=lambda)

    # Number available during each rep (pass)
    for(i in 1:nSites) {
        N[i,] <- rbinom(nReps, M[i], phi)
        }

    # Number observed
    for(i in 1:nSites) {
        for(t in 1:nReps) {
            cp <- c(pA * (1 - pB), pB * (1 - pA), pA * pB)
            cp[4] <- 1 - sum(cp)
            y[i,,t] <- c(rmultinom(1, N[i,t], cp)[1:3])
            }
        }
    return(matrix(y, nSites))
}

str(sim.doub())

# Fit the model

set.seed(4)
y.sim <- sim.doub()
T <- ncol(y.sim) / 3
observer <- matrix(c("A", "B"), 200, T*2, byrow=TRUE)
umf <- unmarkedFrameGMM(y = y.sim,
    obsCovs = list(observer=observer),
    numPrimary=2, type="double")
summary(umf)

m4 <- gmultmix(~1, ~1, ~observer, umf, mixture="NB")
m4

checkEqualsNumeric(coef(m4), c(-0.06998556, 0.77150482, 1.31340048,
    -0.94099309, -1.14215950), tol=1e-5)

backTransform(m4, type="lambda")  # Average abundance per site
backTransform(m4, type="phi")     # Availability
backTransform(linearComb(m4, c(1,0), type="det"))     # obsA detection prob
backTransform(linearComb(m4, c(1,1), type="det"))     # obsB detection prob
backTransform(m4, type="alpha")   # Over-dispersion

# Total pop size
coef(backTransform(m4, type="lambda")) * nrow(y.sim)


pb4 <- parboot(m4, nsim=5, report=1)



nsim <- 500
simout <- matrix(NA, nsim, 5)
colnames(simout) <- c("lambda", "phi", "pA", "pB", "alpha")
for(i in 1:nsim) {
    cat("sim", i, "\n"); flush.console()
    y.sim <- sim.doub()
    T <- ncol(y.sim)/3
    observer <- matrix(c("A", "B"), nrow(y.sim), T*2, byrow=TRUE)
    umf <- unmarkedFrameGMM(y = y.sim, obsCovs=list(observer=observer),
        type="double", numPrimary=T)
    m.sim4 <- gmultmix(~1, ~1, ~observer, umf, mixture="NB")
    e <- coef(m.sim4)
    simout[i,] <- c(exp(e[1]), plogis(e[2:3]), plogis(sum(e[3:4])), exp(e[5]))
    }

hist(simout[,1]); abline(v=1, col=4)
hist(simout[,2]); abline(v=0.6, col=4)
hist(simout[,3]); abline(v=0.8, col=4)
hist(simout[,4]); abline(v=0.6, col=4)
hist(simout[,5]); abline(v=0.5, col=4)












# ------------------------- dependent double observer ------------------



sim.dep.double <- function(nSites=200, numPrimary=2, lambda=1, phi=0.6,
                           pA=0.8, pB=0.6, obsmat, alpha=0.5)
{
    if(numPrimary==1 & phi<1) {
        phi <- 1
        warning("phi has been set to 1 because it can't be estimated when numPrimary=1")
    }

    N <- matrix(NA, nSites, numPrimary)
    y <- array(NA, c(nSites, 2, numPrimary))

    # Abundance at each site
    M <- rnbinom(nSites, size=alpha, mu=lambda)

    # Number available during each rep
    for(i in 1:nSites) {
        N[i,] <- rbinom(numPrimary, M[i], phi)
        }

    if(!all(names(table(obsmat)) == c("A", "B")))
        stop("This function assumes that 'obsmat' is a matrix of A's and B's")
    obsmat <- array(obsmat, c(nSites, 2, numPrimary))

    # Number observed
    for(i in 1:nSites) {
        for(t in 1:numPrimary) {
            if(obsmat[i,1,t]=="A") {
                cp <- c(pA, pB * (1 - pA))
                cp[3] <- 1 - sum(cp)
            }
            if(obsmat[i,1,t]=="B") {
                cp <- c(pB, pA * (1 - pB))
                cp[3] <- 1 - sum(cp)
            }
            y[i,,t] <- c(rmultinom(1, N[i,t], cp)[1:2])
            }
        }
    return(matrix(y, nSites))
}



# piFun

depDoubPiFun <- function(p) {
    M <- nrow(p)
    pi <- matrix(NA, M, 2)
    pi[,1] <- p[,1]
    pi[,2] <- p[,2]*(1-p[,1])
    return(pi)
}

obsToY <- matrix(1, 2, 2)
numPrimary <- 2
obsToY <- kronecker(diag(numPrimary), obsToY)



# Fit the model

set.seed(4)
nSites <- 200
T <- 1
observer <- matrix(c("A", "B", "B", "A"), nSites, T*2, byrow=TRUE)
y.sim <- sim.dep.double(nSites=nSites, numPrimary=T, lambda=3,
                        pA=0.8, pB=0.6, obsmat=observer)
obsToY <- matrix(1, 2, 2)
obsToY <- kronecker(diag(T), obsToY)
umf <- unmarkedFrameGMM(y = y.sim,
    obsCovs = list(observer=observer),
    numPrimary=T, obsToY=obsToY, piFun="depDoubPiFun")
summary(umf)

m5 <- gmultmix(~1, ~1, ~observer-1, umf, mixture="NB")
m5
m6 <- gmultmix(~1, ~1, ~1, umf, mixture="NB")
m6

(pAhat <- plogis(coef(m5, type="det")[1]))
(pBhat <- plogis(coef(m5, type="det")[2]))
1 - (1-pAhat)*(1-pBhat)

plogis(coef(m6, type="det"))

# export data for DOBSERV
dataout <- data.frame(obs1=ifelse(observer[,1]=="A", 1, 2),
                      Species="MALL", y.sim)
head(dataout)
write.csv(dataout, "C:/R/dobserv/simdata.csv", row.names=FALSE)




nsim <- 50
simout <- matrix(NA, nsim, 4)
colnames(simout) <- c("lambda", "phi", "pA", "pB")
for(i in 1:nsim) {
    cat("sim", i, "\n")
    T <- 5
    observer <- matrix(c("A", "B", "B", "A"), nSites, T*2, byrow=TRUE)
    y.sim <- sim.dep.double(nSites=200, alpha=1000, numPrimary=T,
                            obsmat=observer)
    obsToY <- matrix(1, 2, 2)
    obsToY <- kronecker(diag(T), obsToY)
    umf <- unmarkedFrameGMM(y = y.sim, obsCovs=list(observer=observer),
        numPrimary=T, obsToY=obsToY, piFun="depDoubPiFun")
    m.sim5 <- gmultmix(~1, ~1, ~observer-1, umf, mixture="P", se=FALSE)
    e <- coef(m.sim5)
    simout[i,] <- c(exp(e[1]), plogis(e[2:4]))
    cat("   mle =", simout[i,], "\n")
    }

par(mfrow=c(2,2), mai=c(0.8,0.8,0.2,0.2))
hist(simout[,1]); abline(v=1, col=4, lwd=2)
hist(simout[,2]); abline(v=0.6, col=4, lwd=2)
hist(simout[,3]); abline(v=0.8, col=4, lwd=2)
hist(simout[,4]); abline(v=0.6, col=4, lwd=2)

