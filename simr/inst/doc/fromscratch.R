## ---- message=FALSE, warning=FALSE---------------------------------------
library(simr)

## ----options, echo=FALSE, message=FALSE----------------------------------
simrOptions(nsim=100, progress=FALSE)

## ------------------------------------------------------------------------
x <- rep(1:10)
g <- c('a', 'b', 'c')

X <- expand.grid(x=x, g=g)

## ------------------------------------------------------------------------
b <- c(2, -0.1) # fixed intercept and slope
V1 <- 0.5 # random intercept variance
V2 <- matrix(c(0.5,0.05,0.05,0.1), 2) # random intercept and slope variance-covariance matrix
s <- 1 # residual variance

## ------------------------------------------------------------------------
model1 <- makeLmer(y ~ x + (1|g), fixef=b, VarCorr=V1, sigma=s, data=X)
print(model1)
model2 <- makeGlmer(z ~ x + (x|g), family="poisson", fixef=b, VarCorr=V2, data=X)
print(model2)

## ------------------------------------------------------------------------
powerSim(model1, nsim=20)
powerSim(model2, nsim=20)

