library(HyperbolicDist)

### Test gigRawMom
Theta <- c(-0.5,5,2.5)
gigMean(Theta)
gigRawMom(1, Theta)
momIntegrated("gig", order = 1, param = Theta, about = 0)
gigRawMom(2, Theta)
momIntegrated("gig", order = 2, param = Theta, about = 0)
gigRawMom(10, Theta)
momIntegrated("gig", order = 10, param = Theta, about = 0)
gigRawMom(2.5, Theta)
momIntegrated("gig", order = 2.5, param = Theta, about = 0, absolute = TRUE)

### Test gigMom
Theta <- c(-0.5,5,2.5)
m1 <- gigRawMom(1, Theta)
m1
gigMom(1, Theta)
gigMom(2, Theta, m1)
m2 <- momIntegrated("gig", order = 2, param = Theta, about = m1)
m2
gigMom(1, Theta, m1)
gigVar(Theta)
gigMom(3, Theta, m1)
momIntegrated("gig", order = 3, param = Theta, about = m1)
#gigMom(2.5, Theta, m1)

### Test gammaRawMom
shape <- 2
rate <- 3
Theta <- c(shape, rate)
gammaRawMom(1, shape, rate)
momIntegrated("gamma", order = 1, param = Theta, about = 0)
gammaRawMom(2, shape, rate)
momIntegrated("gamma", order = 2, param = Theta, about = 0)
gammaRawMom(10, shape, rate)
momIntegrated("gamma", order = 10, param = Theta, about = 0)

### Test gigMom for gamma case
shape <- 2
rate <- 3
Theta <- c(shape, rate)
m1 <- gammaRawMom(1, shape, rate)
m1
gigMom(1, c(shape, 0, 2*rate))
gigMom(2, c(shape, 0, 2*rate), about = m1)
m2 <- momIntegrated("gamma", order = 2, param = Theta, about = m1)
m2
gigMom(1, c(shape, 0, 2*rate), about = m1)
gigMom(3, c(shape, 0, 2*rate), about = m1)
momIntegrated("gamma", order = 3, param = Theta, about = m1)

### Test gigRawMom for gamma case
library(actuar)
Theta <- c(0.5,0,2.5)
gigRawMom(2, Theta)
mgamma(2, Theta[1], Theta[3]/2)
momIntegrated("gamma", order = 2,
              param = c(Theta[1],Theta[3]/2), about = 0)
### Infinite moment
gigRawMom(-2, Theta)
mgamma(-2, Theta[1], Theta[3]/2)
#momIntegrated("gamma", order = -2,
#              param = c(Theta[1],Theta[3]/2), about = 0)
### Test gigRawMom for inverse gamma case
Theta <- c(-0.5,5,0)
### Infinite moment
gigRawMom(2, Theta)
minvgamma(2, -Theta[1], Theta[2]/2)
#momIntegrated("invgamma", order = 2,
#              param = c(-Theta[1],Theta[2]/2), about = 0)
gigRawMom(-2, Theta)
minvgamma(-2, -Theta[1], Theta[2]/2)
momIntegrated("invgamma", order = -2,
              param = c(-Theta[1],Theta[2]/2), about = 0)

### Test gigMom for gamma case
Theta <- c(0.5,0,2.5)
m1 <- gigMom(1, Theta)
m1
mgamma(1, Theta[1], Theta[3]/2)
momIntegrated("gamma", order = 1,
              param = c(Theta[1],Theta[3]/2), about = 0)
gigMom(2, Theta, m1)
momIntegrated("gamma", order = 2,
              param = c(Theta[1],Theta[3]/2), about = m1)
### Infinite moments
gigMom(-2, Theta)
#momIntegrated("gamma", order = -2,
#              param = c(Theta[1],Theta[3]/2), about = 0)
#gigMom(-2, Theta, m1)
#momIntegrated("gamma", order = -2,
#              param = c(Theta[1],Theta[3]/2), about = m1)

### Test gigRawMom for inverse gamma case
Theta <- c(-0.5,5,0)
### Infinite moments
m1 <- gigMom(1, Theta)
m1
minvgamma(1, -Theta[1], Theta[2]/2)
#momIntegrated("invgamma", order = 1,
#              param = c(-Theta[1],Theta[2]/2), about = 0)
gigMom(2, Theta, m1)
#momIntegrated("invgamma", order = 2,
#              param = c(-Theta[1],Theta[2]/2), about = m1)

gigMom(-2, Theta)
momIntegrated("invgamma", order = -2,
              param = c(-Theta[1],Theta[2]/2), about = 0)
#gigMom(-2, Theta, m1)
momIntegrated("invgamma", order = -2,
              param = c(-Theta[1],Theta[2]/2), about = m1)


