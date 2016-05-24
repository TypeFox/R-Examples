
library(unmarked)
library(rbenchmark)


set.seed(34537)
nSites <- 100
nReps <- 5
covariates <- data.frame(veght=rnorm(nSites),
    habitat=factor(c(rep('A', 50), rep('B', 50))))
psipars <- c(-1, 1, -1)
ppars <- c(1, -1, 0)
X <- model.matrix(~veght+habitat, covariates) # design matrix
psi <- plogis(X %*% psipars)
p <- plogis(X %*% ppars)
y <- matrix(NA, nSites, nReps)
Z <- rbinom(nSites, 1, psi)       # true occupancy state
for(i in 1:nSites) {
    y[i,] <- rbinom(nReps, 1, Z[i]*p[i])
    }
y[1,5,90] <- NA
covariates[2,] <- NA
umf <- unmarkedFrameOccu(y = y, siteCovs = covariates)
head(umf)
summary(umf)

fmR <- occu(~veght ~veght+habitat, umf, engine="C")
fmC <- occu(~veght ~veght+habitat, umf, engine="R")

all.equal(fmR, fmC)
all.equal(coef(fmR), coef(fmC))
all.equal(vcov(fmR), vcov(fmC))


benchmark(occu(~veght ~veght+habitat, umf, engine="C"),
          occu(~veght ~veght+habitat, umf, engine="R"),
          columns=c("test", "elapsed", "relative"),
          replications=500)


