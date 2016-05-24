library(HyperbolicDist)
testVG <- FALSE

### Test ghyp
m1 <- momIntegrated("ghyp", order = 1, param = c(1/2,3,1,1,0), about = 0)
m1
ghypMean(c(1/2,3,1,1,0))
momIntegrated("ghyp", order = 1, param = c(1/2,3,1,1,0), about = m1)
m2 <- momIntegrated("ghyp", order = 2, param = c(1/2,3,1,1,0), about = 0)
m2
m2 - m1^2
momIntegrated("ghyp", order = 2, param = c(1/2,3,1,1,0), about = m1)
momIntegrated("generalized hyperbolic", order = 2, param = c(1/2,3,1,1,0),
              about = m1)
ghypVar(c(1/2,3,1,1,0))


### Test hyperb
m1 <- momIntegrated("hyperb", order = 1, param = c(2,1,1,0), about = 0)
m1
hyperbMean(c(2,1,1,0))
momIntegrated("hyperb", order = 1, param = c(2,1,1,0), about = m1)
m2 <- momIntegrated("hyperb", order = 2, param = c(2,1,1,0), about = 0)
m2
m2 - m1^2
momIntegrated("hyperb", order = 2, param = c(2,1,1,0), about = m1)
momIntegrated("hyperbolic", order = 2, param = c(2,1,1,0),
              about = m1)
hyperbVar(c(2,1,1,0))

### Test gig
m1 <- momIntegrated("gig", order = 1, param = c(1,2,3), about = 0)
m1
gigMean(c(1,2,3))
gigRawMom(1, c(1,2,3))
momIntegrated("gig", order = 1, param = c(1,2,3), about = m1)

m2 <- momIntegrated("gig", order = 2, param = c(1,2,3), about = 0)
m2
m2 - m1^2
momIntegrated("gig", order = 2, param = c(1,2,3), about = m1)
momIntegrated("generalized inverse Gaussian", order = 2, param = c(1,2,3),
              about = m1)
gigVar(c(1,2,3))
gigMom(2, c(1,2,3), about = m1)

### Test gamma
m1 <- momIntegrated("gamma", order = 1, param = c(2,3), about = 0)
m1
shape <- 2
rate <- 3
gigMom(1, c(shape, 0, 2*rate))
momIntegrated("gamma", order = 1, param = c(2,3), about = m1)
m2 <- momIntegrated("gamma", order = 2, param = c(2,3), about = 0)
m2
m2 - m1^2
momIntegrated("gamma", order = 2, param = c(2,3), about = m1)
gigMom(2, c(shape, 0, 2*rate), about = m1)

### Test inverse gamma
m1 <- momIntegrated("invgamma", order = 1, param = c(3,4), about = 0)
m1
shape <- 3
rate <- 1/4
gigMom(order = -1, c(shape, 0, 2*rate))
momIntegrated("invgamma", order = 1, param = c(3,4), about = m1)
m2 <- momIntegrated("invgamma", order = 2, param = c(3,4), about = 0)
m2
gigMom(order = -2, c(shape, 0, 2*rate))
m2 - m1^2
momIntegrated("invgamma", order = 2, param = c(3,4), about = m1)

### Test vg if available
if (testVG){
    library(VarianceGamma)
    m1 <- momIntegrated("vg", order = 1, param = c(2,0.5,0,3), about = 0)
    m1
    vgMean(param = c(2,0.5,0,3))
    momIntegrated("vg", order = 1, param = c(2,0.5,0,3), about = m1)
    m2 <- momIntegrated("vg", order = 2, param = c(2,0.5,0,3), about = 0)
    m2
    m2 - m1^2
    momIntegrated("vg", order = 2, param = c(2,0.5,0,3), about = m1)
    momIntegrated("Variance Gamma", order = 2, param = c(2,0.5,0,3),
                  about = m1)
    vgVar(param = c(2,0.5,0,3))

    m3 <- momIntegrated("vg", order = 3, param = c(2,0.5,0,3), about = 0)
    cm3 <- m3 - 3*m2*m1 + 2*m1^3
    s <- sqrt(m2 - m1^2)
    cm3/(s^3)
    vgSkew(param = c(2,0.5,0,3))
}
