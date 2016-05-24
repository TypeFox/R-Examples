# Copyright 2012-2014 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav

# This file is part of SharpeR.
#
# SharpeR is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SharpeR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SharpeR.  If not, see <http://www.gnu.org/licenses/>.

# env var:
# nb: 
# see also:
# todo:
# changelog: 
#
# Created: 2012.05.19
# Copyright: Steven E. Pav, 2012-2013
# Author: Steven E. Pav
# Comments: Steven E. Pav

# for package roxygen2 help see also
# http://stackoverflow.com/questions/7678284/roxygen2-error-titlerequires-a-value?rq=1

#' Inference on Sharpe ratio and Markowitz portfolio.
#' 
#' @section Sharpe Ratio: 
#'
#' Suppose \eqn{x_i}{xi} are \eqn{n} independent draws of a normal random
#' variable with mean \eqn{\mu}{mu} and variance \eqn{\sigma^2}{sigma^2}.
#' Let \eqn{\bar{x}}{xbar} be the sample mean, and \eqn{s} be
#' the sample standard deviation (using Bessel's correction). Let \eqn{c_0}{c0}
#' be the 'risk free' or 'disastrous rate' of return.  Then
#' \deqn{z = \frac{\bar{x} - c_0}{s}}{z = (xbar - c0)/s} 
#' is the (sample) Sharpe ratio.
#' 
#' The units of \eqn{z} are \eqn{\mbox{time}^{-1/2}}{per root time}.
#' Typically the Sharpe ratio is \emph{annualized} by multiplying by
#' \eqn{\sqrt{d}}{sqrt(d)}, where \eqn{d} is the number of observations
#' per year (or whatever the target annualization epoch.) It is \emph{not}
#' common practice to include units when quoting Sharpe ratio, though
#' doing so could avoid confusion.
#'
#' The Sharpe ratio follows a rescaled non-central t distribution. That
#' is, \eqn{z/K} follows a non-central t-distribution
#' with \eqn{m} degrees of freedom and non-centrality parameter
#' \eqn{\zeta / K}, for some \eqn{K}, \eqn{m} and \eqn{\zeta}{zeta}.
#'
#' We can generalize Sharpe's model to APT, wherein we write
#' \deqn{x_i = \alpha + \sum_j \beta_j F_{j,i} + \epsilon_i,}{
#' x_i = alpha + sum_j beta_j F_j,i + epsilon_i,}
#' where the \eqn{F_{j,i}} are observed 'factor returns', and
#' the variance of the noise term is \eqn{\sigma^2}{sigma^2}.
#' Via linear regression, one can compute estimates \eqn{\hat{\alpha}}{alpha},
#' and \eqn{\hat{\sigma}}{sigma}, and then let the 'Sharpe ratio' be
#' \deqn{z = \frac{\hat{\alpha} - c_0}{\hat{\sigma}}.}{z = (alpha - c0)/sigma.}
#' As above, this Sharpe ratio follows a rescaled t-distribution
#' under normality, \emph{etc.}
#'
#' The parameters are encoded as follows:
#' \itemize{
#' \item \code{df} stands for the degrees of freedom, typically \eqn{n-1}, but
#' \eqn{n-J-1} in general.
#' \item \eqn{\zeta}{zeta} is denoted by \code{zeta}.
#' \item \eqn{d} is denoted by \code{ope}. ('Observations Per Year')
#' \item For the APT form of Sharpe, \code{K} stands for the
#' rescaling parameter.
#' }
#'
#' @section Optimal Sharpe Ratio: 
#'
#' Suppose \eqn{x_i}{xi} are \eqn{n} independent draws of a \eqn{q}-variate
#' normal random variable with mean \eqn{\mu}{mu} and covariance matrix
#' \eqn{\Sigma}{Sigma}. Let \eqn{\bar{x}}{xbar} be the (vector) sample mean, and 
#' \eqn{S} be the sample covariance matrix (using Bessel's correction). Let
#' \deqn{Z(w) = \frac{w^{\top}\bar{x} - c_0}{\sqrt{w^{\top}S w}}}{Z(w) = (w'xbar - c0)/sqrt(w'Sw)}
#' be the (sample) Sharpe ratio of the portfolio \eqn{w}, subject to 
#' risk free rate \eqn{c_0}{c0}.
#'
#' Let \eqn{w_*}{w*} be the solution to the portfolio optimization problem:
#' \deqn{\max_{w: 0 < w^{\top}S w \le R^2} Z(w),}{max {Z(w) | 0 < w'Sw <= R^2},}
#' with maximum value \eqn{z_* = Z\left(w_*\right)}{z* = Z(w*)}.
#' Then 
#' \deqn{w_* = R \frac{S^{-1}\bar{x}}{\sqrt{\bar{x}^{\top}S^{-1}\bar{x}}}}{%
#' w* = R S^-1 xbar / sqrt(xbar' S^-1 xbar)}
#' and
#' \deqn{z_* = \sqrt{\bar{x}^{\top} S^{-1} \bar{x}} - \frac{c_0}{R}}{%
#' z* = sqrt(xbar' S^-1 xbar) - c0/R}
#'
#' The variable \eqn{z_*}{z*} follows an \emph{Optimal Sharpe ratio}
#' distribution. For convenience, we may assume that the sample statistic
#' has been annualized in the same manner as the Sharpe ratio, that is 
#' by multiplying by \eqn{d}, the number of observations per
#' epoch.
#' 
#' The Optimal Sharpe Ratio distribution is parametrized by the number 
#' of assets, \eqn{q}, the number of independent observations, \eqn{n}, the 
#' noncentrality parameter, 
#' \deqn{\zeta_* = \sqrt{\mu^{\top}\Sigma^{-1}\mu},}{zeta* = sqrt(mu' Sigma^-1 mu),}
#' the 'drag' term, \eqn{c_0/R}{c0/R}, and the annualization factor, \eqn{d}.
#' The drag term makes this a location family of distributions, and 
#' by default we assume it is zero.
#' 
#' The parameters are encoded as follows:
#' \itemize{
#' \item \eqn{q} is denoted by \code{df1}.
#' \item \eqn{n} is denoted by \code{df2}.
#' \item \eqn{\zeta_*}{zeta*} is denoted by \code{zeta.s}.
#' \item \eqn{d} is denoted by \code{ope}.
#' \item \eqn{c_0/R} is denoted by \code{drag}.
#' }
#'
#' @section Spanning and Hedging:
#'
#' As above, let 
#' \deqn{Z(w) = \frac{w^{\top}\bar{x} - c_0}{\sqrt{w^{\top}S w}}}{Z(w) = (w'xbar - c0)/sqrt(w'Sw)}
#' be the (sample) Sharpe ratio of the portfolio \eqn{w}, subject to 
#' risk free rate \eqn{c_0}{c0}.
#'
#' Let \eqn{G} be a \eqn{g \times q}{g x q} matrix of 'hedge constraints'. 
#' Let \eqn{w_*}{w*} be the solution to the portfolio optimization problem:
#' \deqn{\max_{w: 0 < w^{\top}S w \le R^2,\,G S w = 0} Z(w),}{max {Z(w) | 0 < w'Sw <= R^2, G S w = 0},}
#' with maximum value \eqn{z_* = Z\left(w_*\right)}{z* = Z(w*)}.
#' Then \eqn{z_*^2}{z*^2} can be expressed as the difference of two squared
#' optimal Sharpe ratio random variables. A monotonic transform takes this
#' difference to the LRT statistic for portfolio spanning, first described by
#' Rao, and refined by Giri. 
#'
#' @section Legal Mumbo Jumbo:
#'
#' SharpeR is distributed in the hope that it will be useful,
#' but WITHOUT ANY WARRANTY; without even the implied warranty of
#' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#' GNU Lesser General Public License for more details.
#'
#' @template etc
#' @template ref-Sharpe
#' @template ref-JW
#' @template ref-Lo
#' @template ref-Opdyke
#' @template ref-LW
#' @template ref-NG1
#'
#' @references
#'
#' Kan, Raymond and Smith, Daniel R. "The Distribution of the Sample Minimum-Variance Frontier."
#' Journal of Management Science 54, no. 7 (2008): 1364--1380.
#' \url{http://mansci.journal.informs.org/cgi/doi/10.1287/mnsc.1070.0852}
#'
#' Kan, Raymond and Zhou, GuoFu. "Tests of Mean-Variance Spanning."
#' Annals of Economics and Finance 13, no. 1 (2012)
#' \url{http://www.aeconf.net/Articles/May2012/aef130105.pdf}
#'
#' Britten-Jones, Mark. "The Sampling Error in Estimates of Mean-Variance 
#' Efficient Portfolio Weights." The Journal of Finance 54, no. 2 (1999):
#' 655--671. \url{http://www.jstor.org/stable/2697722}
#'
#' Silvapulle, Mervyn. J. "A Hotelling's T2-type statistic for testing against 
#' one-sided hypotheses." Journal of Multivariate Analysis 55, no. 2 (1995):
#' 312--319. \url{http://dx.doi.org/10.1006/jmva.1995.1081}
#'
#' Bodnar, Taras and Okhrin, Yarema. "On the Product of Inverse Wishart
#' and Normal Distributions with Applications to Discriminant Analysis 
#' and Portfolio Theory." Scandinavian Journal of Statistics 38, no. 2 (2011):
#' 311--331. \url{http://dx.doi.org/10.1111/j.1467-9469.2011.00729.x}
#'
#' @name SharpeR-package
#' @rdname SharpeR
#' @docType package
#' @title statistics concerning Sharpe ratio and Markowitz portfolio
#' @keywords package
#' @family sr sropt
#' @note The following are still in the works:
#' \enumerate{
#' \item Corrections for standard error based on skew, kurtosis and
#' autocorrelation.
#' \item Tests on Sharpe under positivity constraint. (\emph{c.f.} Silvapulle)
#' \item Portfolio spanning tests.
#' \item Tests on portfolio weights.
#' }
#' This package is maintained as a hobby. 
#'
#' @import matrixcalc 
#' @importFrom stats complete.cases confint cov deviance df dt lm na.omit optimize pchisq pf power.t.test printCoefmat pt qf qnorm qt rchisq rf rnorm rt sd time uniroot vcov
#' @importFrom utils capture.output
NULL

#' @title News for package 'SharpeR':
#'
#' @description 
#'
#' News for package 'SharpeR'
#'
#' \newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
#' \newcommand{\SharpeR}{\CRANpkg{SharpeR}}
#'
#' @section Changes in \SharpeR{} Version 1.1.0 (2016-03-14) :
#' \itemize{
#' \item fix sr_vcov on array input.
#' \item add SRIC method.
#' \item add SRIC to print.sropt.
#' \item change predint output to matrix.
#' }
#'
#' @section Changes in \SharpeR{} Version 1.0.0 (2015-06-18) :
#'
#' \itemize{
#' \item sane version numbers.
#' \item unpaired k sample test of Sharpe.
#' \item rely on same for unpaired 2 sample test.
#' \item prediction intervals for Sharpe based on upsilon.
#' \item more tests.
#' }
#'
#' @section Changes in \SharpeR{} Version 0.1501 (2014-12-06) :
#'
#' \itemize{
#' \item fix inference of mark frequency from e.g. xts objects.
#' \item add rlambdap.
#' }
#'
#' @section Changes in \SharpeR{} Version 0.1401 (2014-01-05) :
#'
#' \itemize{
#' \item fix second moment asymptotic covariance.
#' \item add confidence distribution functions for sr, sr.opt.
#' }
#'
#' @section Changes in \SharpeR{} Version 0.1310 (2013-10-30) :
#'
#' \itemize{
#' \item inverse second moment asymptotic covariance.
#' }
#'
#' @section Changes in \SharpeR{} Version 0.1309 (2013-09-20) :
#'
#' \itemize{
#' \item spanning/hedging tests.
#' \item sr equality test via callback variance covariance computation.
#' \item split vignette in two.
#' }
#' 
#' @section Changes in \SharpeR{} Version 0.1307 (2013-05-30) :
#'
#' \itemize{
#' \item proper d.f. in sr objects with different nan fill.
#' \item restore vignette.
#' }
#'
#' @section \SharpeR{} Initial Version 0.1306 (2013-05-21) :
#' \itemize{
#' \item put on CRAN
#' }
#'
#' @name SharpeR-NEWS
#' @rdname NEWS
NULL

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
