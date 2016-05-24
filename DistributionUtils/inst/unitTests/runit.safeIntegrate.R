### Unit tests of function safeIntegrate

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make

test.safeIntegrate <- function()
{
  ## Purpose: Level 1 test of safeIntegrate
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Local, Date: 03 Feb 2010, 16:25

  checkEquals(safeIntegrate(dnorm, -1.96, 1.96)$value,
              0.9500042097036, tolerance = 1.0e-11)
  checkEquals(safeIntegrate(dnorm, 1.96, 1.96)$value,
              0, tolerance = 1.0e-11)
  checkEquals(safeIntegrate(dnorm, -Inf, Inf)$value,
              0.9999999997942, tolerance = 1.0e-11)
  checkEquals(safeIntegrate(dnorm, -Inf, -Inf)$value,
              0, tolerance = 1.0e-11)
  checkEquals(safeIntegrate(dnorm, Inf, Inf)$value,
              0, tolerance = 1.0e-11)
  return()
}

