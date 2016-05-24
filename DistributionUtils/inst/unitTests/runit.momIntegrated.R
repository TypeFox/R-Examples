### Unit tests of function momIntegrated

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make

test.momIntegrated <- function()
{
  ## Gamma Distribution
  ## Raw moments of gamma
  rawMom <- numeric(8)
  gammaMom <- function(order, shape, scale){
    gMom <- (scale^order)*gamma(shape + order)/gamma(shape)
    return(gMom)
  }

  ## Calculate moments for particular gamma
  shape <- 2
  rate <- 3
  scale <- 1
  rawMom <- sapply(1:8, gammaMom, shape = shape, scale = scale)
  ## Central moments, gamma
  centralMom <- momChangeAbout("all", rawMom, 0, rawMom[1])
  ## Moments about new value
  new <- 1
  newMom <- momChangeAbout("all", rawMom, 0, new)

  ## Check integrated moments from gamma
  ## Raw moments
  m1 <- momIntegrated("gamma", order = 1,
                      shape = shape, scale = scale, about = 0)
  m8 <- momIntegrated("gamma", order = 8,
                      shape = shape, scale = scale, about = 0)
  checkEquals(rawMom[1], m1)
  checkEquals(rawMom[8], m8)

  ## Central moments
  cm1 <- momIntegrated("gamma", order = 1,
                       shape = shape, scale = scale, about = m1)
  cm8 <- momIntegrated("gamma", order = 8,
                       shape = shape, scale = scale, about = m1)
  checkEquals(centralMom[1], cm1)
  checkEquals(centralMom[8], cm8)

  ## Moments about new
  nm1 <- momIntegrated("gamma", order = 1,
                       shape = shape, scale = scale,
                       about = new)
  nm8 <- momIntegrated("gamma", order = 8,
                       shape = shape, scale = scale,
                       about = new)
  checkEquals(newMom[1], nm1)
  checkEquals(newMom[8], nm8)

  ## Normal Distribution
  ## Central moment for normal distribution
  centralMom <- numeric(8)
  normalMom <- function(order, mean, sd){
    if (order%%2 == 0){
      nMom <- sd^order * sqrt(2^order/pi) * gamma((order - 1)/2 + 1)
    }else{
      nMom <- 0
    }
    return(nMom)
  }

  ## Calculate moments for particular normal distribution
  mean <- 2
  sd <- 1
  centralMom <- sapply(1:8, normalMom, mean = mean, sd = sd)
  ## Raw moments, normal distribution
  rawMom <- momChangeAbout("all", centralMom, mean, 0)
  ## Moments about new value
  new <- 1
  newMom <- momChangeAbout("all", centralMom, mean, new)

  ## Check integrated moments from normal
  ## Raw moments
  m1 <- momIntegrated("norm", order = 1,
                      mean = mean, sd = sd, about = 0)
  m8 <- momIntegrated("norm", order = 8,
                      mean = mean, sd = sd, about = 0)
  checkEquals(rawMom[1], m1)
  checkEquals(rawMom[8], m8)

  ## Central moments
  cm1 <- momIntegrated("norm", order = 1,
                       mean = mean, sd = sd, about = m1)
  cm8 <- momIntegrated("norm", order = 8,
                       mean = mean, sd = sd, about = m1)
  checkEquals(centralMom[1], cm1)
  checkEquals(centralMom[8], cm8)

  ## Moments about new
  nm1 <- momIntegrated("norm", order = 1,
                       mean = mean, sd = sd, about = new)
  nm8 <- momIntegrated("norm", order = 8,
                       mean = mean, sd = sd, about = new)
  checkEquals(newMom[1], nm1)
  checkEquals(newMom[8], nm8)

  return()
}


