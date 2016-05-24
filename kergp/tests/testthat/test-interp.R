library(testthat)
library(kergp)


## copied from DiceKriging
branin <- function (x) {
    x1 <- x[1] * 15 - 5
    x2 <- x[2] * 15
    (x2 - 5/(4 * pi^2) * (x1^2) + 5/pi * x1 - 6)^2 +
        10 * (1 - 1/(8 * pi)) * cos(x1) + 10
}

d <- 2; n <- 16
design.fact <- expand.grid(x1 = seq(0, 1, length = 4),
                           x2 = seq(0, 1, length = 4))
y <- apply(design.fact, 1, branin) 

## kriging model 1 : 2D Gaussian kernel
myGaussFunVec <- function(X1, X2, par) { 
    ## X1, X2 : matrices with same number of column d (dimension)
    n <- nrow(X1)
    d <- ncol(X1)     
    SS2 <- 0
    
    for (j in 1:d){
        Aj <- outer(X1[, j], X2[, j], "-")
        Hj2 <- (Aj / par[1])^2
        SS2 <- SS2 + Hj2
    }
    D2 <- exp(-SS2)
    kern <- par[2]*D2
    D1 <- 2 * kern * SS2 / par[1] 
    attr(kern, "gradient") <- list(theta = D1,  sigma2 = D2)
    
    return(kern)
}

myGaussVec <- covMan(
    kernel = myGaussFunVec,
    hasGrad = TRUE,
    acceptMatrix = TRUE,
    d = 2,
    parLower = c(theta = 0.0, sigma2 = 0.0),
    parUpper = c(theta = +Inf, sigma2 = +Inf),
    parNames = c("theta", "sigma2"),
    label = "my Gaussian kernel"
)

myCoef <- c(theta = 0.5, sigma2 = 10000)

m <- gp(formula = y ~ 1, data = data.frame(y = y, design.fact),
        inputs = names(design.fact), cov = myGaussVec,
        parCovIni = myCoef, 
        beta = 0, varNoiseLower = 0, varNoiseUpper = 0,
        parCovLower = myCoef, parCovUpper = myCoef) 

mNoise <- gp(formula = y ~ 1, data = data.frame(y = y, design.fact),
              inputs = names(design.fact), cov = myGaussVec,
              parCovIni = myCoef, 
              beta = 0, varNoiseLower = 1000, varNoiseUpper = 1000,
              parCovLower = myCoef, parCovUpper = myCoef) 

p <- predict(m, newdata = design.fact, type = "UK", forceInterp = TRUE)
pNoise <- predict(mNoise, newdata = design.fact, type = "UK", forceInterp = TRUE)

precision <- 1e-10
## the following tests should work with it,
## since the computations are analytical

test_that(desc="Kriging mean (no noise), on the design points", 
          expect_true(max(abs(p$mean - y)) < precision))

test_that(desc="Kriging mean with noise, on the design points", 
          expect_true(max(abs(pNoise$mean - y)) < precision))


