
library(unmarked)
library(rbenchmark)


set.seed(34537)
nSites <- 100
nReps <- 5
covariates <- data.frame(veght=rnorm(nSites),
    habitat=factor(c(rep('A', 50), rep('B', 50))))
lampars <- c(-1, 1, -1)
ppars <- c(1, -1, 0)
X <- model.matrix(~veght+habitat, covariates) # design matrix
lam <- exp(X %*% lampars)
p <- plogis(X %*% ppars)
y <- matrix(NA, nSites, nReps)
N <- rnbinom(nSites, size=1, mu=lam)       # true abund
for(i in 1:nSites) {
    y[i,] <- rbinom(nReps, N[i], p[i])
    }
y[1,5,90] <- NA
covariates[2,] <- NA
umf <- unmarkedFramePCount(y = y, siteCovs = covariates)
head(umf)
summary(umf)

fm1R <- pcount(~veght ~veght+habitat, umf, engine="R", K=50)
fm1C <- pcount(~veght ~veght+habitat, umf, engine="C", K=50)

all.equal(fm1R, fm1C)
all.equal(coef(fm1R), coef(fm1C))
all.equal(vcov(fm1R), vcov(fm1C))


benchmark(pcount(~veght ~veght+habitat, umf, engine="C", K=30),
          pcount(~veght ~veght+habitat, umf, engine="R", K=30),
          columns=c("test", "elapsed", "relative"),
          replications=100)


fm2R <- pcount(~veght ~veght+habitat, umf, engine="R", K=50, mixture="NB")
fm2C <- pcount(~veght ~veght+habitat, umf, engine="C", K=50, mixture="NB")

all.equal(fm2R, fm2C)
all.equal(coef(fm2R), coef(fm2C))
all.equal(vcov(fm2R), vcov(fm2C))

