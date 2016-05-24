### Unit tests of function besselRatio

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make

test.besselRatio <- function()
{
  ## Purpose: Level 1 test of besselRatio
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Local, Date: 25 Jan 2010, 14:17

  nus <- c(0:5, 10, 20)
  x <- seq(1, 4, length.out = 11)
  k <- 3

  ## check calculation of ratios agrees with calculation via besselK
  bK <- matrix(nrow = length(nus), ncol = length(x))
  bRatio <- matrix(nrow = length(nus), ncol = length(x))
  compare <- matrix(nrow = length(nus), ncol = length(x))
  unitMatrix <- matrix(1, nrow = length(nus), ncol = length(x))
  
  for (i in 1:length(nus)){
    for (j in 1:length(x)) {
      bK[i,j] <- besselK(x[j], nus[i] + k)/besselK(x[j], nus[i])
      bRatio[i,j] <- besselRatio(x[j], nus[i],
                                 orderDiff = k)

      compare[i,j] <- bK[i,j]/bRatio[i,j]
    }
  }

  checkEquals(compare, unitMatrix)  

  ## check exponential scaling works properly
  raw <- matrix(nrow = length(nus), ncol = length(x))
  scaled <- matrix(nrow = length(nus), ncol = length(x))
  compare <- matrix(nrow = length(nus), ncol = length(x))
  unitMatrix <- matrix(1, nrow = length(nus), ncol = length(x))

  for (i in 1:length(nus)){
    for (j in 1:length(x)) {
      raw[i,j] <- besselRatio(x[j], nus[i],
                              orderDiff = k)
      scaled[i,j] <- besselRatio(x[j], nus[i],
                                 orderDiff = k, useExpScaled = 1)
      compare[i,j] <- raw[i,j]/scaled[i,j]
    }
  }

  checkEquals(compare, unitMatrix) 

  return()
}

