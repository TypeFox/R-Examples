## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


### Extreme-value copulas ######################################################

evCopula <- function(family, param = NA_real_, dim = 2L, ...) {
  familiesImplemented <- c("galambos", "gumbel", "huslerReiss", "tawn", "tev")
  fam <- pmatch(family, familiesImplemented, -1)
  if (fam == -1)
    stop("Valid family names are ",
         paste(familiesImplemented, collapse=", "))

  switch(fam,
	 galambosCopula	  (param),
	 gumbelCopula	  (param),
	 huslerReissCopula(param),
	 tawnCopula	  (param),
	 tevCopula	  (param),
	 stop("family ", fam, "not yet available, at least via evCopula()"))
}

tailIndexEvCopula <- function(copula) {
  lower <- 0
  upper <- 2 - 2 * A(copula, 0.5)
  c(lower=lower, upper=upper)
}


hdensity <- function(copula, z) {
  A <- A(copula, z)
  ders <- dAdu(copula, z)
  1 + (1 - 2 * z) * ders$der1 / A + z * (1 - z) * (ders$der2 * A - ders$der1^2) / A^2
}

revCopula <- function(n, copula) {
  ## Caperaa, Fougeres, and Genest (2000, Journal of Multivariate Analysis)
  ## This algorithm has low efficiency for high dependence, in which case,
  ## hdensity is pretty much centered around 0.5.
  ## In particular, it generates peculiar numbers for galambosCopula with
  ## alpha >= 30. Don't how to solve yet.

  M <- hdensity(copula, 0.5) ## maximum obtained at 0.5 for symmetric copula
  z <- rep(NA, n)
  ndone <- 0
  while (TRUE) {
    ucand <- runif(n)
    accept <- runif(n) <= hdensity(copula, ucand) / M
    accept <- accept & (!is.na(accept)) ## hdensity can be NA at some ucand
    ngood <- sum(accept)
    if (ngood == 0) next
    ngood <- min(ngood, n - ndone)
    z[ndone + 1:ngood] <- ucand[accept][1:ngood]
    ndone <- ndone + ngood
    if (ndone == n) break
  }
  ders <- dAdu(copula, z)
  Az <- A(copula, z)
  pz <- z * (1 - z) * ders$der2 / hdensity(copula, z) / Az
  w <- numeric(n)
  mix1 <- runif(n) <= pz
  nmix1 <- sum(mix1)
  if (any( mix1)) w[ mix1] <- runif(nmix1)
  if (any(!mix1)) w[!mix1] <- runif(n - nmix1) * runif(n - nmix1)
  ## CFG (2000, p.39)
  l.A <- log(w)/Az
  exp(cbind(z * l.A, (1 - z) * l.A))
}

### These one-dimensional numerical integrations are quite accurate.
### They are much better than two-dimensional integration based on function adapt.

tauEvCopula <- function(copula) {
  integrand <- function(x) x * (1 - x) / A(copula, x) * dAdu(copula, x)$der2
  integrate(integrand, 0, 1)$value
}

rhoEvCopula <- function(copula) {
  integrand <- function(x) 1 / (A(copula, x) + 1)^2
  12 * integrate(integrand, 0, 1)$value - 3
}

setMethod("tailIndex", signature("evCopula"), tailIndexEvCopula)
setMethod("rCopula", signature("numeric", "evCopula"), revCopula)
setMethod("tau", signature("evCopula"), tauEvCopula)
setMethod("rho", signature("evCopula"), rhoEvCopula)

