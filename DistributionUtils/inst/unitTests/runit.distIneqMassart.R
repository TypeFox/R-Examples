### Unit tests of function distIneqMassart

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make

test.distIneqMassart <- function()
{
  ## Purpose: Level 1 test of distIneqMassart
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date:  9 Feb 2010, 12:11

  ## Normal (default)
  checkTrue(distIneqMassart()$check)
  ## Specify parameter values of normal
  checkTrue(distIneqMassart(mean = 1, sd = 2)$check)
  
  ## Gamma distribution has no default value for shape
  checkTrue(distIneqMassart("gamma", shape = 1)$check)

  return()
}
