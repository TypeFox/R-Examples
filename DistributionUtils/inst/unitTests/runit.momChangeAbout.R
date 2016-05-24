### Unit tests of function momChangeAbout

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make

test.momChangeAbout <- function()
{
  ## Purpose: Level 1 test of momChangeAbout
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date:  4 Feb 2010, 14:34


  ## Gamma distribution
  k <- 4
  shape <- 2
  scale <- 1
  old <- 0
  new <- shape*scale         # central moments
  sampSize <- 10000
  x <- rgamma(sampSize, shape, scale = scale)


  ## Sample moments
  s4new <- mean((x - new)^k)
  s3new <- mean((x - new)^3)


  ## Calculate 1st to 4th raw moments
  m <- numeric(k)
  for (i in 1:k){
    m[i] <- gamma(shape + i)/gamma(shape)
  }


  ## Calculate 4th moment about new
  m4new <- momChangeAbout(k, m, old, new)
  ## Calculate 3rd about new
  m3new <- momChangeAbout(3, m, old, new)

  ## Calculate standard errors for gamma
  rawMom <- numeric(8)
  gammaMom <- function(order, shape, scale){
    gMom <- (scale^order)*gamma(shape + order)/gamma(shape)
    return(gMom)
  }
  rawMom <- sapply(1:8, gammaMom, shape = shape, scale = scale)
  centralMom <- momChangeAbout("all", rawMom, 0, rawMom[1])
  s4SE <- momSE(4, sampSize, centralMom)
  s3SE <- momSE(3, sampSize, centralMom[1:6])
  ## Compare with sample values
  s4tol <- qnorm(0.995)*s4SE
  s3tol <- qnorm(0.995)*s3SE
  checkTrue(abs(s4new - m4new) < s4tol)
  checkTrue(abs(s3new - s3new) < s3tol)

  return()
}
