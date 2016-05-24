### Unit tests of function integrateDens

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make

test.integrateDens <- function()
{
  ## Purpose: Level 1 test of integrateDens
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date: 08 Feb 2010, 21:07

  ## Normal (default)
  totPr <- integrateDens()$value
  checkEquals(totPr, 1)
  totPr <- integrateDens("norm", mean = 3, sd = 2)$value
  checkEquals(totPr, 1)
  half <- integrateDens("norm", lower = 3, mean = 3, sd = 2)$value
  checkEquals(half, 1/2)

  ## t
  totPr <- integrateDens("t", df = 4)$value
  checkEquals(totPr, 1)
  half <- integrateDens("t", lower = 0, df = 4)$value
  checkEquals(half, 1/2)

  ## Weibull
  totPr <- integrateDens("weibull", shape = 1)
  checkEquals(half, 1/2)

  return()
}


