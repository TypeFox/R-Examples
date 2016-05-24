### Unit tests of functions  in sampleMoments.R

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make

test.sampleMoments <- function()
{
  ## Purpose: Level 1 test of sampleMoments
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date: 03 Feb 2010, 23:00

  ## Check sample moments from normal
  sampSize <- 1000
  sigma <- 2
  x <- rnorm(sampSize, sd = sigma)
  distSkew <- 0
  distKurt <- 0
  sampSkew <- skewness(x)
  sampKurt <- kurtosis(x)
  diffSkew <- abs(distSkew - sampSkew)
  diffKurt <- abs(distKurt - sampKurt)

  ## Calculate tolerances
  ## Central moments, standard normal
  mom <- c(0,sigma^2,0,3*sigma^4,0,15*sigma^6,0,105*sigma^8)
  s3SE <- momSE(3, sampSize, mom[1:6])
  s4SE <- momSE(4, sampSize, mom)
  tolSkew <- qnorm(0.995)*s3SE/(sigma^3)
  tolKurt <- qnorm(0.995)*s4SE/(sigma^4)
  checkTrue(diffSkew < tolSkew)
  checkTrue(diffKurt < tolKurt)

  ## Check sample moments from gamma
  sampSize <- 10000
  shape <- 8
  scale <- 2
  set.seed(123)
  x <- rgamma(sampSize, shape = shape, scale = scale)
  sampM3 <- skewness(x)*sd(x)^3
  sampM4 <- (kurtosis(x) + 3)*sd(x)^4
  ## Calculate tolerances
  ## Raw moments of gamma
  rawMom <- numeric(8)
  gammaMom <- function(order, shape, scale){
    gMom <- (scale^order)*gamma(shape + order)/gamma(shape)
    return(gMom)
  }
  rawMom <- sapply(1:8, gammaMom, shape = shape, scale = scale)
  ## Central moments, gamma
  centralMom <- momChangeAbout("all", rawMom, 0, rawMom[1])
  distSD <- centralMom[2]
  diffM3 <- abs(centralMom[3] - sampM3)
  diffM4 <- abs(centralMom[4] - sampM4)
  s3SE <- momSE(3, sampSize, centralMom[1:6])
  s4SE <- momSE(4, sampSize, centralMom)
  ## gamma is very skewed: need to allow large tolerance
  tolM3 <- qnorm(0.999)*s3SE
  tolM4 <- qnorm(0.995)*s4SE
  checkTrue(diffM3 < tolM3,
            msg = paste("tolM3 = ", tolM3, " diffM3 = ", diffM3))
  checkTrue(diffM4 < tolM4,
            msg = paste("tolM4 = ", tolM4, " diffM4 = ", diffM4))

  return()
}

graphicstest.sampleMoments <- function()
{
  ## Purpose: Check sample moment functions graphically
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date:  4 Feb 2010, 12:07

  ## Standard normal
  sampSizes <- c(5,10,20,50,100,200,500,1000,2000,5000,10000)
  results <- matrix(nrow = length(sampSizes), ncol = 4,
                    dimnames = list(NULL,c("tsk","ssk","tku","sku")))
  results <- as.data.frame(results)
  results[, 1] <- 0
  results[, 3] <- 0
  for (i in (1:length(sampSizes))){
    x <- rnorm(sampSizes[i])
    results[i, 2] <- skewness(x)
    results[i, 4] <- kurtosis(x)
  }

  ## plot results
  with(results, {plot(sampSizes, ssk,
                     xlab = "Sample Size", ylab = "Sample Skewness")
                 lines(sampSizes, tsk)
                 })
  title("Sample and Theoretical Skewness: Standard Normal")
  with(results, {plot(sampSizes, sku,
                     xlab = "Sample Size", ylab = "Sample Kurtosis")
                 lines(sampSizes, tku)})
  title("Sample and Theoretical Kurtosis: Standard Normal")


  ## Gamma(1,1)
  sampSizes <- c(5,10,20,50,100,200,500,1000,2000,5000,10000)
  results <- matrix(nrow = length(sampSizes), ncol = 4,
                    dimnames = list(NULL,c("tsk","ssk","tku","sku")))
  results <- as.data.frame(results)
  results[, 1] <- 2
  results[, 3] <- 6
  for (i in (1:length(sampSizes))){
    x <- rgamma(sampSizes[i], shape = 1)
    results[i, 2] <- skewness(x)
    results[i, 4] <- kurtosis(x)
  }

  ## plot results
  graphicsOutput <- paste(pathReport, "sampleMoments.pdf", sep = "")
  cat("Graphics output in file ", graphicsOutput, "\n")
  pdf(file = graphicsOutput)
  print(pathReport)
  with(results, {plot(sampSizes, ssk,
                     xlab = "Sample Size", ylab = "Sample Skewness")
                 lines(sampSizes, tsk)
                 })
  title("Sample and Theoretical Skewness: Gamma(1,1)")
  with(results, {plot(sampSizes, sku,
                     xlab = "Sample Size", ylab = "Sample Kurtosis")
                 lines(sampSizes, tku)})
  title("Sample and Theoretical Kurtosis: Gamma(1,1)")

  dev.off()

  return()
}


