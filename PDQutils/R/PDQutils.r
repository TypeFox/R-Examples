# Copyright 2015-2015 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav
#
# This file is part of PDQutils.
#
# PDQutils is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PDQutils is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with PDQutils.  If not, see <http://www.gnu.org/licenses/>.

# Created: 2015.01.31
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav
# Comments: Steven E. Pav

#' PDQ Functions via Gram-Charlier, Edgeworth, and Cornish Fisher Approximations
#' 
#' @section Gram Charlier and Edgeworth Expansions:
#'
#' Given the raw moments of a probability distribution, we can approximate the probability 
#' density function, or the cumulative distribution function, via a Gram-Charlier 'A'
#' expansion on the standardized distribution.
#'
#' Suppose \eqn{f(x)}{f(x)} is the probability density of some random
#' variable, and let \eqn{F(x)}{F(x)} be the cumulative distribution function.
#' Let \eqn{He_j(x)}{He_j(x)} be the \eqn{j}{j}th probabilist's Hermite
#' polynomial. These polynomials form an orthogonal basis, with respect to the
#' function \eqn{w(x)}{w(x)} of the Hilbert space of functions which are square
#' \eqn{w}{w}-weighted integrable. The weighting function is 
#' \eqn{w(x) = e^{-x^2/2} = \sqrt{2\pi}\phi(x)}{w(x) = e^{-x^2/2} = sqrt(2pi) phi(x)}.
#' The orthogonality relationship is
#' \deqn{\int_{-\infty}^{\infty} He_i(x) He_j(x) w(x) \mathrm{d}x = \sqrt{2\pi} j! \delta_{ij}.}{integral_-inf^inf He_i(x) He_j(x) w(x) dx = sqrt(2pi)j!dirac_ij.}
#'
#' Expanding the density \eqn{f(x)}{f(x)} in terms of these polynomials in the
#' usual way (abusing orthogonality) one has
#' \deqn{f(x) \approx \sum_{0\le j} \frac{He_j(x)}{j!} \phi(x) \int_{-\infty}^{\infty} f(z) He_j(z) \mathrm{d}z.}{f(x) = sum_{0 <= j} (He_j(x)/j!) phi(x) integral_-inf^inf f(z) He_j(z) dz.}
#' The cumulative distribution function is 'simply' the integral of this
#' expansion. Abusing certain facts regarding the PDF and CDF of the normal
#' distribution and the probabilist's Hermite polynomials, the CDF has
#' the representation
#' \deqn{F(x) = \Phi(x) - \sum_{1\le j} \frac{He_{j-1}(x)}{j!} \phi(x) \int_{-\infty}^{\infty} f(z) He_j(z) \mathrm{d}z.}{F(x) = Phi(x) - sum_{1 <= j} (He_{j-1}(x)/j!) phi(x) integral_-inf^inf f(z) He_j(z) dz.}
#'
#' These series contain coefficients defined by the probability distribution 
#' under consideration. They take the form
#' \deqn{c_j = \frac{1}{j!}\int_{-\infty}^{\infty} f(z) He_j(z) \mathrm{d}z.}{c_j = (1/j!) integral_-inf^inf f(z) He_j(z) dz.}
#' Using linearity of the integral, these coefficients are easily computed in
#' terms of the coefficients of the Hermite polynomials and the raw, uncentered
#' moments of the probability distribution under consideration. Note that it may be the
#' case that the computation of these coefficients suffers from bad numerical
#' cancellation for some distributions, and that an alternative formulation
#' may be more numerically robust.
#'
#' @section Generalized Gram Charlier Expansions:
#'
#' The Gram Charlier 'A' expansion is most appropriate for random variables
#' which are vaguely like the normal distribution. For those which are like
#' another distribution, the same general approach can be pursued. One needs
#' only define a weighting function, \eqn{w}, which is the density of the
#' 'parent' probability distribution, then find polynomials,
#' \eqn{p_n(x)}{p_n(x)} which are orthogonal with respect to \eqn{w} over
#' its support. One has
#' \deqn{f(x) = \sum_{0\le j} p_j(x) w(x) \frac{1}{h_j} \int_{-\infty}^{\infty} f(z) p_j(z) \mathrm{d}z.}{f(x) = sum_{0 <= j} p_j w(x) (1/h_j) integral_-inf^inf f(z) p_j(z) dz.}
#' Here \eqn{h_j}{h_j} is the normalizing constant:
#' \deqn{h_j = \int w(z)p_j^2(z)\mathrm{d}z.}{h_j = integral w(z)p_j^2(z) dz.}
#' One must then use facts about the orthogonal polynomials to approximate
#' the CDF.
#'
#' Another approach to arrive at the same computation is described by 
#' Berberan-Santos.
#'
#' @section Cornish Fisher Approximation:
#'
#' The Cornish Fisher approximation is the Legendre
#' inversion of the Edgeworth expansion of a distribution, but ordered
#' in a way that is convenient when used on the mean of a number of
#' independent draws of a random variable. 
#'
#' Suppose \eqn{x_1, x_2, \ldots, x_n}{x_1, x_2, ..., x_n} are \eqn{n} independent 
#' draws from some probability distribution. 
#' Letting 
#' \deqn{X = \frac{1}{\sqrt{n}} \sum_{1 \le i \le n} x_i,}{X = (x_1 + x_2 + ... x_n) / sqrt(n),}
#' the Central Limit Theorem assures us that, assuming finite variance, 
#' \deqn{X \rightarrow \mathcal{N}(\sqrt{n}\mu, \sigma),}{X ~~ N(sqrt(n) mu, sigma),}
#' with convergence in \eqn{n}.
#'
#' The Cornish Fisher approximation gives a more detailed picture of the
#' quantiles of \eqn{X}{X}, one that is arranged in decreasing powers of
#' \eqn{\sqrt{n}}{sqrt(n)}. The quantile function is the function \eqn{q(p)}{q(p)} 
#' such that \eqn{P\left(X \le q(p)\right) = q(p)}{P(x <= q(p)) = p}. The
#' Cornish Fisher expansion is 
#' \deqn{q(p) = \sqrt{n}\mu + \sigma \left(z + \sum_{3 \le j} c_j f_j(z)\right),}{q(p) = sqrt{n}mu + sigma (z + sum_{3 <= j} c_j f_j(z)),}
#' where \eqn{z = \Phi^{-1}(p)}{z = qnorm(p)}, and \eqn{c_j}{c_j} involves
#' standardized cumulants of the distribution of \eqn{x_i}{x_i} of order
#' \eqn{j} and higher. Moreover, the \eqn{c_j}{c_j} feature decreasing powers
#' of \eqn{\sqrt{n}}{sqrt(n)}, giving some justification for truncation.
#' When \eqn{n=1}{n=1}, however, the ordering is somewhat arbitrary.
#'
#' @section Legal Mumbo Jumbo:
#'
#' PDQutils is distributed in the hope that it will be useful,
#' but WITHOUT ANY WARRANTY; without even the implied warranty of
#' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#' GNU Lesser General Public License for more details.
#'
#' @template etc
#'
#' @template ref-AS269
#' @template ref-Jaschke
#' @template ref-Blinnikov
#' @template ref-bsantos
#'
#' @import orthopolynom moments
#' @importFrom stats coef dbeta dgamma dnorm pgamma pnorm qnorm runif
#'
#' @name PDQutils
#' @rdname PDQutils
#' @docType package
#' @title PDQ Functions via Gram Charlier, Edgeworth, and Cornish Fisher Approximations
#' @keywords package
#' 
#' @note 
#'
#' This package is maintained as a hobby. 
#'
NULL

#' @title News for package \sQuote{PDQutils}:
#'
#' @description
#'
#' News for package \sQuote{PDQutils}.
#'
#' \newcommand{\pkg}{#1}
#' \newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
#' \newcommand{\PDQutils}{\CRANpkg{PDQutils}}
#'
#' @section \PDQutils{} Version 0.1.3 (2016-01-04) :
#' \itemize{
#' \item Package maintenance--no new features.
#' }
#'
#' @section \PDQutils{} Version 0.1.2 (2015-06-15) :
#' \itemize{
#' \item Generalized Gram Charlier expansions.
#' \item bugfixes.
#' }
#'
#' @section \PDQutils{} Version 0.1.1 (2015-02-26) :
#' \itemize{
#' \item Edgeworth expansions.
#' }
#'
#' @section \PDQutils{} Initial Version 0.1.0 (2015-02-14) :
#' \itemize{
#' \item first CRAN release.
#' }
#'
#' @name PDQutils-NEWS
#' @rdname NEWS
NULL

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
