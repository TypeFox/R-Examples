### Testing nlMean
test.nlMean <- function() {
  param <- c(0, 1, 2, 3)
  testMean <- nlMean(param = param)

  ddist <- function(x, order, param, about) {
    (x - about)^order*dnl(x, param = param)
  }

  momMean <- integrate(ddist, -30, 30, param = param, order = 1,
                       about = 0, subdivisions = 1000,
                       rel.tol = .Machine$double.eps^0.5)[[1]]

  checkEquals(testMean, momMean)
}


### Testing nlVar
test.nlVar <- function() {
  param <- c(0, 1, 2, 3)
  testVar <- nlVar(param = param)
  mn <- nlMean(param = param)

  ddist <- function(x, order, param, about) {
    (x - about)^order*dnl(x, param = param)
  }

  momVar <- integrate(ddist, -30, 30, param = param, order = 2,
                      about = mn, subdivisions = 1000,
                      rel.tol = .Machine$double.eps^0.5)[[1]]

  checkEquals(testVar, momVar)
}


### Testing nlSkew
test.nlSkew <- function() {
  param <- c(0, 1, 2, 3)
  testSkew <- nlSkew(param = param)
  mn <- nlMean(param = param)

  ddist <- function(x, order, param, about) {
    (x - about)^order*dnl(x, param = param)
  }

  m3 <- integrate(ddist, -30, 30, param = param, order = 3,
                  about = mn, subdivisions = 1000,
                  rel.tol = .Machine$double.eps^0.5)[[1]]
  sigma3 <- nlVar(param = param)^(3/2)
  momSkew <- m3/sigma3

  checkEquals(testSkew, momSkew)
}


### Testing nlKurt
test.nlKurt <- function() {
  param <- c(0, 1, 2, 3)
  testKurt <- nlKurt(param = param)
  mn <- nlMean(param = param)

  ddist <- function(x, order, param, about) {
    (x - about)^order*dnl(x, param = param)
  }

  m4 <- integrate(ddist, -30, 30, param = param, order = 4,
                  about = mn, subdivisions = 1000,
                  rel.tol = .Machine$double.eps^0.5)[[1]]
  sigma4 <- nlVar(param = param)^2
  momKurt <- m4/sigma4 - 3

  checkEquals(testKurt, momKurt)
}
