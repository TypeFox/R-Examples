##
## runit-normal.R - Unit test for normal type DF
##
## Author:
##  Olaf Mersmann <olafm@statistik.tu-dortmund.de>
##

d1min <- 0.022750131948179209
d2min <- 0.150831468618899983
dmax <- 0.977249868051820791

test.normMin <- function() {
  d <- normMin(-1, 1)
  checkEqualsNumeric(d(-1), dmax)
  checkEqualsNumeric(d(0), 0.5)
  checkEqualsNumeric(d(1), d1min)
}

test.normMinArgs <- function() {
  checkException(normMin(2, 1))
  checkException(normMax(1, 1))
  checkException(normMax(2, 1))
  checkException(normTarget(2, 1, 0))
  checkException(normTarget(1, 2, 2))
  checkException(normTarget(1, 1, 2))
}

test.normMax <- function() {
  d <- normMax(-1, 1)
  checkEqualsNumeric(d(-1), d1min)
  checkEqualsNumeric(d(0), 0.5)
  checkEqualsNumeric(d(1), dmax)
}

test.normTarget <- function() {
  d <- normTarget(-1, 0, 1)
  checkEqualsNumeric(d(-1), d2min)
  checkEqualsNumeric(d(0), dmax)
  checkEqualsNumeric(d(1), d2min)

  z <- seq(-2, 0, by=.2)
  checkEqualsNumeric(d(z), d(-z))  
}

checkNormDesireDensity <- function(d) {
  message("")
  s <- seq(0, dmax, by=0.05)
  for (i in 1:10) {
    mean <- runif(1, -2, 2)
    sd <- runif(1, 0, 2)
    message(sprintf("Y ~ N(%f, %f)", mean, sd))
    checkTrue(all(ddesire(s, d, mean, sd) >= 0))
    checkEqualsNumeric(ddesire(-1, d, mean, sd), 0)
    checkEqualsNumeric(ddesire(dmax + .Machine$double.eps, d, mean, sd), 0)
  }
}

test.normMinDensity <- function()
  checkNormDesireDensity(normMin(-1,1))

test.normMaxDensity <- function()
  checkNormDesireDensity(normMax(-1, 1))

test.normTargetDensity <- function()
  checkNormDesireDensity(normTarget(-1, 0, 1))
