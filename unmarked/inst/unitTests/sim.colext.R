
library(unmarked)
library(RUnit)


# ----------------------------- simulate ---------------------------------


sim <- function(nSites=100, nReps=5, nYears=5, psi=0.5, gamma=0.2,
                epsilon=0.8, p=0.4)
{

    y <- array(NA, c(nSites, nReps, nYears))
    Z <- matrix(NA, nSites, nYears)

    phi <- 1-epsilon

    Z[,1] <- rbinom(nSites, 1, psi)
    for(t in 2:nYears) {
        muZ <- Z[,t-1] * phi + (1 - Z[,t-1]) * gamma
        Z[,t] <- rbinom(nSites, 1, muZ)
        }
    for(j in 1:nReps)
        for(t in 1:nYears)
            y[,j,t] <- rbinom(nSites, 1, Z[,t]*p)

    y <- matrix(y, nSites, nReps*nYears)
    return(y)
}

# sim()


# ------------------------------- unmarked -------------------------------


library(unmarked)
library(RUnit)

set.seed(3)
nYears <- 5
sim1 <- sim(nYears=nYears)
umf <- unmarkedMultFrame(y = sim1, numPrimary = nYears)

(m <- colext(~1, ~1, ~1, ~1, umf, control=list(trace=T, REPORT=1)))

backTransform(m, type="psi")
backTransform(m, type="col")
backTransform(m, type="ext")
backTransform(m, type="det")

checkEqualsNumeric(coef(m),
                   c(-0.1047112, -1.3000613, 1.5203993, -0.3634747),
                   tol=1e-5)









# Covariates

library(unmarked)

nSites <- 100
nReps <- 4
nYears <- 5

set.seed(3454)
x1 <- rnorm(nSites)
x2 <- matrix(rnorm(nSites*nYears), nSites, nYears)

psi <- plogis(-1 + 1*x1)
epsilon <- plogis(-2 + 1*x2)
phi <- 1-epsilon
gamma <- 0.4
p <- 0.3

y <- array(NA, c(nSites, nReps, nYears))
Z <- matrix(NA, nSites, nYears)

Z[,1] <- rbinom(nSites, 1, psi)
for(t in 2:nYears) {
    muZ <- Z[,t-1] * phi[,t-1] + (1 - Z[,t-1]) * gamma
    Z[,t] <- rbinom(nSites, 1, muZ)
}
for(j in 1:nReps)
    for(t in 1:nYears)
    y[,j,t] <- rbinom(nSites, 1, Z[,t]*p)

y <- matrix(y, nSites, nReps*nYears)

#x1[1] <- NA
#x2[1,2] <- NA
#y[1,1:4] <- NA

umf <- unmarkedMultFrame(y=y, siteCovs=data.frame(x1=x1),
                         yearlySiteCovs=list(x2=x2),
                         numPrimary=nYears)

summary(umf)

(m2 <- colext(~x1, ~1, ~x2, ~1, umf))

checkEqualsNumeric(coef(m2), c(-0.7440136, 0.9246523, -0.1577932,
                               -1.9102425,  1.2321186, -0.6352266),
                   tol=1e-6)












# Missing values, no covariates






set.seed(3)
nYears <- 5
sim3 <- sim(nYears=nYears, psi=0.5, gamma=0.2, epsilon=0.8, p=0.4)
sim3[1,1] <- NA
sim3[2,1:5] <- NA
sim3[3,6:10] <- NA
umf3 <- unmarkedMultFrame(y = sim3, numPrimary = nYears)

(m3 <- colext(~1, ~1, ~1, ~1, umf3, control=list(trace=T, REPORT=1)))

backTransform(m3, type="psi")
backTransform(m3, type="col")
backTransform(m3, type="ext")
backTransform(m3, type="det")

checkEqualsNumeric(coef(m),
                   c(-0.1047112, -1.3000613, 1.5203993, -0.3634747),
                   tol=1e-5)























# Missing values and Covariates

library(unmarked)

nSites <- 100
nReps <- 4
nYears <- 5

set.seed(3454)
x1 <- rnorm(nSites)
x2 <- matrix(rnorm(nSites*nYears), nSites, nYears)

psi <- plogis(-1 + 1*x1)
epsilon <- plogis(-2 + 1*x2)
phi <- 1-epsilon
gamma <- 0.4
p <- 0.3

y <- array(NA, c(nSites, nReps, nYears))
Z <- matrix(NA, nSites, nYears)

Z[,1] <- rbinom(nSites, 1, psi)
for(t in 2:nYears) {
    muZ <- Z[,t-1] * phi[,t-1] + (1 - Z[,t-1]) * gamma
    Z[,t] <- rbinom(nSites, 1, muZ)
}
for(j in 1:nReps)
    for(t in 1:nYears)
    y[,j,t] <- rbinom(nSites, 1, Z[,t]*p)

y <- matrix(y, nSites, nReps*nYears)



y[1,1] <- NA
y[2,1:5] <- NA
y[3,6:10] <- NA

x1[c(4,6)] <- NA

x2[7,1] <- NA
x2[8,3] <- NA

umf4 <- unmarkedMultFrame(y=y, siteCovs=data.frame(x1=x1),
                         yearlySiteCovs=list(x2=x2),
                         numPrimary=nYears)

summary(umf4)

(m4 <- colext(~x1, ~1, ~x2, ~1, umf4))
























trace(unmarked:::handleNA, browser, browser, signature="unmarkedMultFrame")
untrace(unmarked:::handleNA, signature="unmarkedMultFrame")

trace(unmarked:::getDesign, browser, browser,
      signature="unmarkedMultFrame")
untrace(unmarked:::getDesign, signature="unmarkedMultFrame")


debugonce(unmarked:::colext.fit)
