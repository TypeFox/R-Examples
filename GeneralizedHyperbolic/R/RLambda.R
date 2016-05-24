### RLambda & SLambda - functions used in calculation of the
### theoretical mean and variance of a hyperbolic distribution
### (as given in Barndorff-Nielsen & Blaesild, 1983).
RLambda <- function(zeta, lambda = 1) {
  besselK(zeta, lambda + 1) / besselK(zeta, lambda)
}

SLambda <- function(zeta, lambda = 1) {
  (besselK(zeta, lambda + 2) * besselK(zeta, lambda) -
   besselK(zeta, lambda + 1)^2) / besselK(zeta, lambda)^2
}

MLambda <- function(zeta, lambda = 1) {
  besselK(zeta, lambda + 2) / besselK(zeta, lambda)
}


### WLambda1, WLambda2, WLambda3 and WLambda4 from
### Barndorff-Nielsen and Blaesild (1981), p. 41
WLambda1 <- function(zeta, lambda = 1) {
  RLambda(zeta, lambda)
}

WLambda2 <- function(zeta, lambda = 1) {
  -RLambda(zeta, lambda)^2 + 2 * (lambda + 1) * RLambda(zeta, lambda) / zeta + 1
}

WLambda3 <- function(zeta, lambda = 1) {
  2 * RLambda(zeta, lambda)^3 - 6 * (lambda + 1) * RLambda(zeta, lambda)^2 /
  zeta + (4 * (lambda + 1) * (lambda + 2) / zeta^2 - 2) * RLambda(zeta, lambda) +
  2 * (lambda + 2) / zeta
}

WLambda4 <- function(zeta, lambda = 1) {
  -6 * RLambda(zeta, lambda)^4 + 24 * (lambda + 1) * RLambda(zeta, lambda)^3 /
  zeta + (-4 * (lambda + 1) * (7 * lambda + 11) / zeta^2 + 8) * RLambda(zeta, lambda)^2 +
  (8 * (lambda + 1) * (lambda + 2) * (lambda + 3) / zeta^3 - 4 * (4 * lambda + 5) / zeta) *
  RLambda(zeta, lambda) + 4 * (lambda + 2) * (lambda + 3) / zeta^2 - 2
}

### Theoretical skewness and kurtosis
### Barndorff-Nielsen and Blaesild (1981), p.42
gammaLambda1 <- function(hyperbPi, zeta, lambda = 1) {
  (hyperbPi^3 * WLambda3(zeta, lambda) + 3 * hyperbPi * WLambda2(zeta, lambda) / zeta) /
  (hyperbPi^2 * WLambda2(zeta, lambda) + WLambda1(zeta) / zeta)^(3 / 2)
}

gammaLambda2 <- function(hyperbPi, zeta, lambda = 1) {
  (hyperbPi^4 * WLambda4(zeta, lambda) +
  6 * hyperbPi^2 * WLambda3(zeta, lambda) / zeta +
  3 * WLambda2(zeta, lambda) / zeta^2) /
  (hyperbPi^2 * WLambda2(zeta, lambda) + WLambda1(zeta) / zeta)^2
}
