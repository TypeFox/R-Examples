
##' Prediction of a Ornstein-Uhlenbeck Gaussian process.
##'
##' The process y(t) is assumed to be given by yprime(t) = alpha *
##' y(t) + epsilon(t) where epsilon(t) is white noise with variance
##' sigma^2. The stationnary distribution is normal with zero mean and
##' variance sigma^2 / 2 / alpha. Due to the Markov property the
##' prediction has a simple closed form.
##'
##' For the sake of simplicity, it is expected that 't' does not embed
##' dupplicates and that no 'newt' element is exactly at an observation
##' time.
##' 
##' @title Prediction of a Ornstein-Uhlenbeck Gaussian process 
##' @param newt new 
##' @param t observation time
##' @param y observation value
##' @param alpha the OU parameter (inverse of a duration). Must be
##' positive. The larger alpha is, the faster the process returns to
##' mu.
##' @param sigma sd of the white noise. When a 'newt' element is
##' far away from the observation time, the prediction boils
##' down to 0 with a prediction variance of sigma^2 / 2 / alpha.
##' @param plot For graphic diag.
##' @return A list with the mean and the variance.
##' @author yves
##' 
predictOUG <- function(newt, t, y,
                       alpha = 0.2, sigma = 1,
                       plot = TRUE) {
    t <- sort(t)
    n <- length(t)
    newn <- length(newt)
    yL <- yR <- rep(0, newn)
    phiL <- phiR <- phi <- rep(0, newn)
    tL <- tR <- t
    
    ## int <- tnew %in% t
    ind <- findInterval(newt, t)

    ## 'tL' and 'tR' are the vector of length 'newn' containing the
    ## observation time at the left and at the right of the t[i]
    hasL <- (ind != 0)
    indL <- ind[hasL]
    tL[hasL] <- t[indL]
    phiL[hasL] <- exp(-alpha * (newt[hasL] - tL[hasL]))
    yL[hasL] <- y[indL] 
    
    hasR <- (ind < n)
    indR <- ind[hasR] + 1
    tR[hasR] <- t[indR]
    phiR[hasR] <- exp(alpha * (newt[hasR] - tR[hasR]))
    yR[hasR] <- y[indR]
    
    delta <- tR - tL
    phi[hasL & hasR] <- exp(-alpha * delta[hasL & hasR])
    ySmooth <- phiL * (yL - phi * yR) + phiR * (yR - phi * yL)
    ySmooth <- ySmooth / ( 1 - phi^2)
    varSmooth <- (1 + phi^2 - phiL^2 - phiR^2) / (1 - phi^2) * sigma^2 / 2 / alpha 
    ## variant: !hasL flags the elts in 'newt' with newt < t[1]
    if (FALSE) {
        varSmooth[!hasL] <- (1 - phiR[!hasL]^2) * sigma^2 / 2 / alpha
        varSmooth[!hasR] <- (1 - phiL[!hasR]^2) * sigma^2 / 2 / alpha
    }
    ## print(data.frame(newt, yL, yR, phiL, phiR, varL, varR))
    if (plot) {
        tLim <- range(t, newt)
        yLim <- range(y, 0) ## use 'mu' for a general rule
        d <- diff(tLim)
        tLim <- tLim + c(-0.05, 0.05) * d
        d <- diff(yLim)
        yLim <- yLim + c(-0.05, 0.05) * d
   
        plot(t, y, type = "n", xlim = tLim, ylim = yLim)
        a <- 2 * sqrt(varSmooth)
        segments(x0 = newt, x1 = newt, y0 = ySmooth - a,
                 y1 = ySmooth + a, col = "lavender")
        abline(h = 0, col = "SpringGreen3", lwd = 2) ## mu
        points(newt, ySmooth, type = "p", pch = 16, col = "SteelBlue2", cex = 1.1)
        abline(v = t, col = "orange")
        points(t, y, type = "p", pch = 21, lwd = 3, col = "orangered", bg = "yellow")
    }
        
    list(mean = ySmooth,
         sd = sqrt(varSmooth),
         var = varSmooth)
    
}

n <- 4
newn <- 100 ## 100000 would not be a problem here
t <- sort(runif(n))
y <- sin(2 * pi * t)
newt <- sort(runif(newn, min = -0.3, max = 2.3))

alpha <- 4; sigma <- 0.2
res <- predictOUG(newt = newt, t = t, y = y, alpha = alpha, sigma = sigma, 
                  plot = FALSE)

## comparison with DiceKriging formulas
mykm <- kmData(y ~ 1, data = data.frame(y, t),
               inputnames = "t",
               covtype = "exp",
               coef.trend = 0,
               coef.cov = 1 / alpha,
               coef.var = sigma^2 / alpha / 2)

pred <- predict(mykm, newdata = data.frame(t = newt), type = "SK")
# points(newt, pred$mean, pch = 15, cex = 0.4, col = "Chartreuse2")

precision <- 1e-10  # the following tests should work with it, since the computations are analytical
test_that(desc="Kriging mean, 'exp' cov. fun.", expect_true(max(abs(res$mean-pred$mean)) < precision))
test_that(desc="Kriging mean, 'exp' cov. fun.", expect_true(max(abs(res$sd-pred$sd)) < precision))


