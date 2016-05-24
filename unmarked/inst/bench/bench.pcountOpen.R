
library(unmarked)
library(rbenchmark)


set.seed(34537)
nSites <- 100
T <- 10
covariates <- data.frame(veght=rnorm(nSites),
    habitat=factor(c(rep('A', 50), rep('B', 50))))
lampars <- c(-1, 0, 0)
gampars <- c(-1, 0, 0)
ompars <- c(2, -1, 0)
ppars <- c(1, 0, 0)
X <- model.matrix(~veght+habitat, covariates) # design matrix
lam <- exp(X %*% lampars)
gam <- exp(X %*% gampars)
om <- plogis(X %*% ompars)
p <- plogis(X %*% ppars)
y <- N <- matrix(NA, nSites, T)
N[,1] <- rnbinom(nSites, size=1, mu=lam)       # true abund
for(t in 2:T) {
    S <- rbinom(nSites, N[,t-1], om)
    G <- rpois(nSites, gam)
    N[,t] <- S+G
    }
y[] <- rbinom(nSites*T, N, p)
#y[1,5] <- NA
#covariates[2,] <- NA
umf <- unmarkedFramePCO(y = y, siteCovs = covariates, numPrimary=T)
head(umf)
summary(umf)

system.time(fm <- pcountOpen(~veght+habitat, ~1, ~veght, ~veght, umf,
                              K=20, control=list(trace=TRUE, REPORT=1),
                              se=FALSE)) # 105 (was 722)!!




benchmark(pcountOpen(~1, ~1, ~veght, ~1, umf, K=15, se=FALSE,
                     starts=c(-1,-1,2,-1,1),
                     control=list(trace=TRUE, REPORT=1)),
          columns=c("test", "elapsed", "relative"),
          replications=5) # 94.25 was 275.43



st2 <- system.time(fm2C <- pcountOpen(~1, ~1, ~veght, ~1, umf, K=15,
                                      control=list(trace=TRUE, REPORT=1),
                                      se=TRUE)) # 51.46 was 155.5


