### Unit tests of function momSE

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make

test.momSE <- function()
{
  ## Purpose: Level 1 test of momSE
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date: 5 Feb 2010, 11:53

  ## Moments of the normal distribution, mean 1, variance 4
  mu <- 1
  sigma <- 2
  n <- 100
  mom <- c(0,sigma^2,0,3*sigma^4,0,15*sigma^6,0,105*sigma^8)
  ## standard error of sample variance
  checkEquals(momSE(2, 100, mom[1:4]), sqrt(2*sigma^4/n))

  ## standard error of sample central third moment
  checkEquals(momSE(3, 100, mom[1:6]), sqrt(6*sigma^6/n))

  ## standard error of sample central fourth moment
  checkEquals(momSE(4, 100, mom), sqrt(96*sigma^8/n))

  ## Moments of the gamma distribution, shape 1, scale 2
  shape <- 1
  scale <- 2
  n <- 1000
  rawMom <- numeric(8)
  gammaMom <- function(order, shape, scale){
    gMom <- (scale^order)*gamma(shape + order)/gamma(shape)
    return(gMom)
  }
  rawMom <- sapply(1:8, gammaMom, shape = shape, scale = scale)
  centralMom <- momChangeAbout("all", rawMom, 0, rawMom[1])

  stErrors <- numeric(4)
  stErrors[1] <- sqrt(centralMom[2]/n)
  for (i in (2:4)){
    stErrors[i] <- momSE(i, n, centralMom[1:(2*i)])
  }
  stErrors
  shouldBe <- c(0.06324555320337,0.35777087639997,
                3.71806401235912,60.10550723519435)


  ## compare with previously calculated answer
  checkEqualsNumeric(stErrors, shouldBe)



  return()
}

