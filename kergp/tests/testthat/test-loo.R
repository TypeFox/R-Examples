library(testthat)
library(kergp)

x <- c(0, 0.4, 0.6, 0.8, 1)
y <- c(-0.3, 0, -0.8, 0.5, 0.9)
n <- length(x)

theta <- 0.01; sigma <- 3;

myGaussFunVec <- function(X1, X2, par) { 
    ## X1, X2 : matrices with same number of column d (dimension)
    n <- nrow(X1)
    d <- ncol(X1)     
    SS2 <- 0
    
    for (j in 1:d){
        Aj <- outer(X1[ , j], X2[ , j], "-")
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
    d = 1,
    parLower = c(theta = 0.0, sigma2 = 0.0),
    parUpper = c(theta = +Inf, sigma2 = +Inf),
    parNames = c("theta", "sigma2"),
    label = "my Gaussian kernel"
)

design <- data.frame(x = x)
m <- gp(formula = y ~ ., data = data.frame(y = y, design),
        inputs = names(design), cov = myGaussVec,
        compGrad = TRUE, checkNames = FALSE,
        parCovIni = c(0.5, 0.5),
        varNoiseLower = 1e-4, varNoiseUpper = 1e-1,
        parCovLower = c(1e-5, 1e-5), parCovUpper = c(Inf, Inf)) 

##
list.type <- list("SK", "SK", "UK", "UK")
list.trend.reestim <- list(FALSE, TRUE, FALSE, TRUE)

precision <- 1e-6

for (b in 1:length(list.type)) {
    if (list.trend.reestim[[b]]) {
        beta <- NULL
    } else {
        beta <- c(-1, 2)
    }
    
    case <- paste(list.type[[b]], "; trend.reestim=", list.trend.reestim[[b]])
    print(case)
    
    loo.mean <- loo.sd <- rep(NA, n)
    for (i in 1:n){
        
        ## downdate of m (should be implemented as a method in a future stage)
        mloo <- m
        mloo$y <- mloo$y[-i]
        mloo$X <- mloo$X[-i, , drop = FALSE]
        mloo$F <- mloo$F[-i, , drop = FALSE]
        mloo$dim$n <- mloo$dim$n - 1
        o <- gls(mloo$covariance, y = mloo$y, X = mloo$X, 
                 F = mloo$F, varNoise = mloo$varNoise,
                 beta = beta, checkNames = FALSE)
        mloo$L <- o$L
        mloo$eStar <- o$eStar
        mloo$sseStar <- o$sseStar
        if (is.null(beta)){  # trend reestimation
            mloo$betaHat <- o$betaHat
            mloo$FStar <- o$FStar
            mloo$RStar <- o$RStar
        } else {  # same trend
            mloo$FStar <- m$FStar[-i, , drop = FALSE]
            mloo$RStar <- qr.R(qr(mloo$FStar))
        }
        ## end downdate of m
        
        p <- predict(mloo, newdata = data.frame(x = x[i]), 
                     type = list.type[[b]], checkNames = FALSE,
                     forceInterp = TRUE)
        loo.mean[i] <- p$mean
        loo.sd[i] <- p$sd
    }
    
    loo <- influence(m, type = list.type[[b]], trend.reestim = list.trend.reestim[[b]])
    
    test_that(desc = paste("Check LOO mean for", list.type[[b]]),
              expect_true(max(abs(loo$mean - loo.mean) / loo.mean) < precision))
    
    test_that(desc = paste("Check LOO sd for", list.type[[b]]),
              expect_true(max(abs(loo$sd - loo.sd) / loo.sd) < precision))

    sum((loo$mean - loo.mean)^2)
    sum((loo$sd - loo.sd)^2)
}
  
