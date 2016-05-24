### Unit tests of the function gigFit
### By David Cusack, 18/10/2010

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make

### Unit test for fitting the Generalized Inverse Gaussian Distribution
test.gigFit <- function() {
  ## Purpose: Level 1 test of gigFit
  data(gigParam)

  errorFit <- 10
  nFit <- 10
  testParam <- matrix(c(1,1,1), nrow = 1)
  #testParam <- gigSmallParam
  n <- 1000

  for (i in 1:NROW(testParam)) {
    param <- testParam[i,]
    #cat(param, "\n", file = "/home/dscott/Packages/GeneralizedHyperbolicDevel/GeneralizedHyperbolic/inst/unitTests/Parameters.txt", append = TRUE)
    chi <- param[1]
    psi <- param[2]
    lambda <- param[3]
    chiHat <- vector(length = nFit)
    psiHat <- vector(length = nFit)
    lambdaHat <- vector(length = nFit)

    for (j in 1:nFit) {
      x <- rgig(n = n, param = param)
      gigFitReturn <- gigFit(x)$param
      chiHat[j] <- gigFitReturn[1]
      psiHat[j] <- gigFitReturn[2]
      lambdaHat[j] <- gigFitReturn[3]
    }
    chiHatAvg <- mean(chiHat)
    psiHatAvg <- mean(psiHat)
    lambdaHatAvg <- mean(lambdaHat)


    checkTrue(abs(chiHatAvg - chi) < errorFit,
              msg = paste(param[1], param[2], param[3], chiHatAvg))
    checkTrue(abs(psiHatAvg - psi) < errorFit,
              msg = paste(param[1], param[2], param[3], psiHatAvg))
    checkTrue(abs(lambdaHatAvg - lambda) < errorFit,
              msg = paste(param[1], param[2], param[3], lambdaHatAvg))
  }
}

### Graphical Test for Generalized Inverse Gaussian Distribution
graphicstest.gigFit <- function() {
  data(gigParam)
  n <- 1000

  smallparam <- gigSmallParam
  smallchi <- smallparam[, 1]
  smallpsi <- smallparam[, 2]
  smalllambda <- smallparam[, 3]
  sden <- 0
  pdf(file = "Histogram and Log Histograms of gigFit.pdf",
      encoding = "ISOLatin1", height = 7, width = 10)
  par(mfrow = c(1, 2), oma = c(5, 5, 5, 5))
  for (i in 1 : nrow(smallparam))
  {
    x <- rgig(n, param = smallparam[i, ])
    fitparam <- try(gigFit(x, startValues = "MoM")$param, silent = TRUE)
    c <- class(fitparam)
    if (c == "try-error") {
      stop(print(smallparam[i, ]))
    }
    sde <- density(x, bw = 0.1)$y
    sden <- dgig(x, param = fitparam)
    hist(x, freq = FALSE, breaks = 20, main = "", xlab = "sample")
    mtext(expression(bold("Graph Test of gigFit")), line = 3.5, cex = 1.15)
    mtext(bquote(paste(chi ==.(smallchi[i]), ",",
                       psi ==.(smallpsi[i]), ",",
                       lambda ==.(smalllambda[i]), sep = "")),
          line = 2.25, cex = 1.15)
    curve(dgig(x, param = fitparam), add = TRUE, col = "red")
    logHist(x, main = "", breaks = 20, htype = "h")
    mtext(expression(bold("Log Graph Test of gigFit")),
          line = 3.5, cex = 1.15)
    mtext(bquote(paste(chi ==.(smallchi[i]), ",",
                       psi ==.(smallpsi[i]), ",",
                       lambda ==.(smalllambda[i]), sep = "")),
          line = 2.25, cex = 1.15)
    curve(log(dgig(x, param = smallparam[i, ])), add = TRUE, col = "red")
    i = i + 1
  }
  dev.off()
}

