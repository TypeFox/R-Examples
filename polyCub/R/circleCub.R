################################################################################
### Part of the R package "polyCub".
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Copyright (C) 2013-2014 Sebastian Meyer
### Time-stamp: <[circleCub.R] by SM Die 06/05/2014 10:02 (CEST)>
###
### Special cases of cubature over circular domains (center, r)
################################################################################


##' Integration of the Isotropic Gaussian Density over Circular Domains
##'
##' This function calculates the integral of the bivariate, isotropic Gaussian
##' density (i.e. \eqn{\Sigma} = \code{sd^2*diag(2)}) over circular domains via
##' the cumulative distribution function of the (non-central) Chi-Squared
##' distribution (\code{pchisq}), cp. Formula 26.3.24 in Abramowitz and Stegun
##' (1972).
##' 
##' @references
##' Abramowitz, M. and Stegun, I. A. (1972).
##' Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical
##' Tables. New York: Dover Publications.
##' @param center numeric vector of length 2 (center of the circle).
##' @param r numeric (radius of the circle). Several radii may be supplied. 
##' @param mean numeric vector of length 2
##'             (mean of the bivariate Gaussian density).
##' @param sd numeric (common standard deviation of the isotropic
##'           Gaussian density in both dimensions).
##' @return The integral value (one for each supplied radius).
##' @note The non-centrality parameter of the evaluated chi-squared distribution
##' equals the squared distance between the \code{mean} and the 
##' \code{center}. If this becomes too large, the result becomes inaccurate, see
##' \code{\link{pchisq}}.
##' @keywords math spatial
##' @importFrom stats pchisq
##' @export
##' @examples
##' circleCub.Gauss(center=c(1,2), r=3, mean=c(4,5), sd=6)
##'
##' if (requireNamespace("mvtnorm") && gpclibPermit()) {
##'   ## compare with cubature over a polygonal approximation of a circle
##'   disc.poly <- spatstat::disc(radius=3, centre=c(1,2), npoly=32)
##'   polyCub.exact.Gauss(disc.poly, mean=c(4,5), Sigma=6^2*diag(2))
##' }

circleCub.Gauss <- function (center, r, mean, sd)
{
    stopifnot(isScalar(sd), length(center) == 2, length(mean) == 2)
    pchisq((r/sd)^2, df=2, ncp=sum(((center-mean)/sd)^2))
}
