# Copyright 2014-2014 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav
#
# This file is part of sadists.
#
# sadists is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# sadists is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with sadists.  If not, see <http://www.gnu.org/licenses/>.

# Created: 2014.02.13
# Copyright: Steven E. Pav, 2014
# Author: Steven E. Pav
# Comments: Steven E. Pav

#' Some Additional Distributions.
#'
#' A collection of distributions which can be approximated via
#' Edgeworth and Cornish-Fisher expansions
#' 
#' @section Sum of (non-central) chi-square to powers:
#'
#' Let \eqn{X_i \sim \chi^2\left(\delta_i, \nu_i\right)}{X_i ~ chi^2(delta_i, v_i)}
#' be independently distributed non-central chi-squares, where \eqn{\nu_i}{v_i}
#' are the degrees of freedom, and \eqn{\delta_i}{delta_i} are the
#' non-centrality parameters.  
#' Let \eqn{w_i} and \eqn{p_i} be given constants. Suppose
#' \deqn{Y = \sum_i w_i X_i^{p_i}.}{Y = sum w_i (X_i)^(p_i).}
#' Then \eqn{Y}{Y} follows a weighted sum of chi-squares to power distribution. 
#' The special case where all the \eqn{p_i}{p_i} are one is a 'sum of
#' chi-squares' distribution; 
#' The special case where all the \eqn{p_i}{p_i} are one half is a 'sum of
#' chis' distribution; 
#'
#' @section Lambda Prime:
#'
#' Introduced by Lecoutre, the lambda prime distribution
#' finds use in inference on the Sharpe ratio under normal
#' returns.
#' Suppose \eqn{y \sim \chi^2\left(\nu\right)}{y ~ x^2(v)}, and
#' \eqn{Z}{Z} is a standard normal. 
#' \deqn{T = Z + t \sqrt{y/\nu}}{T = Z + t sqrt(y/v)}
#' takes a lambda prime distribution with parameters 
#' \eqn{\nu, t}{v, t}.
#' A lambda prime random variable can be viewed as a confidence
#' variable on a non-central t because 
#' \deqn{t = \frac{Z' + T}{\sqrt{y/\nu}}}{t = (Z' + T)/sqrt(y/v)}
#'
#'
#' @section Upsilon:
#'
#' The upsilon distribution generalizes the lambda prime to the
#' case of the sum of multiple chi variables. That is,
#' suppose \eqn{y_i \sim \chi^2\left(\nu_i\right)}{y_i ~ x^2(v_i)}
#' independently and independently of \eqn{Z}{Z}, a standard normal. 
#' Then 
#' \deqn{T = Z + \sum_i t_i \sqrt{y_i/\nu_i}}{T = Z + sum_i t_i sqrt(y_i/v_i)}
#' takes an upsilon distribution with parameter vectors
#' \eqn{[\nu_1, \nu_2, \ldots, \nu_k]', [t_1, t_2, ..., t_k]'}{<v_1, v_2, ..., v_k>, <t_1, t_2, ..., t_k>}.
#'
#' The upsilon distribution is used in certain tests of
#' the Sharpe ratio for independent observations.
#'
#' @section K Prime:
#'
#' Introduced by Lecoutre, the K prime family of distributions generalize
#' the (singly) non-central t, and lambda prime distributions. 
#' Suppose \eqn{y \sim \chi^2\left(\nu_1\right)}{y ~ x^2(v1)}, and
#' \eqn{x \sim t \left(\nu_2, a\sqrt{y/\nu_1}/b\right)}{x ~ t(v2,(a/b) sqrt(y/v1))}.
#' Then the random variable
#' \deqn{T = b x}{T = b x}
#' takes a K prime distribution with parameters 
#' \eqn{\nu_1, \nu_2, a, b}{v1, v2, a, b}. In Lecoutre's terminology,
#' \eqn{T \sim K'_{\nu_1, \nu_2}\left(a, b\right)}{T ~ K'_v1,v2(a,b)}
#'
#' Equivalently, we can think of
#' \deqn{T = \frac{b Z + a \sqrt{\chi^2_{\nu_1} / \nu_1}}{\sqrt{\chi^2_{\nu_2} / \nu_2}}}{T = (bZ + a sqrt(chi2_v1/v1)) / sqrt(chi2_v2/v2)}
#' where \eqn{Z} is a standard normal, and the normal and the (central) chi-squares are
#' independent of each other. When \eqn{a=0}{a=0} we recover
#' a central t distribution; 
#' when \eqn{\nu_1=\infty}{v1=inf} we recover a rescaled non-central t distribution;
#' when \eqn{b=0}{b=0}, we get a rescaled square root of a central F
#' distribution; when \eqn{\nu_2=\infty}{v2=inf}, we recover a 
#' Lambda prime distribution.
#'
#' @section Doubly Noncentral t:
#'
#' The doubly noncentral t distribution generalizes the (singly)
#' noncentral t distribution to the case where the numerator is
#' the square root of a scaled noncentral chi-square distribution.
#' That is, if 
#' \eqn{X \sim \mathcal{N}\left(\mu,1\right)}{X ~ N(u,1)} independently
#' of \eqn{Y \sim \chi^2\left(k,\theta\right)}{Y ~ x^2(k,theta)}, then
#' the random variable
#' \deqn{T = \frac{X}{\sqrt{Y/k}}}{T = X / sqrt(Y/k)}
#' takes a doubly non-central t distribution with parameters
#' \eqn{k, \mu, \theta}{k, mu, theta}.
#'
#' @section Doubly Noncentral F:
#'
#' The doubly noncentral F distribution generalizes the (singly)
#' noncentral F distribution to the case where the numerator is
#' a scaled noncentral chi-square distribution.
#' That is, if 
#' \eqn{X \sim \chi^2\left(n_1,\theta_1\right)}{X ~ x^2(n1,theta1)} independently 
#' of \eqn{Y \sim \chi^2\left(n_2,\theta_2\right)}{Y ~ x^2(n2,theta2)}, then
#' the random variable
#' \deqn{T = \frac{X / n_1}{Y / n_2}}{T = (X/n1) / (Y/n2)}
#' takes a doubly non-central F distribution with parameters
#' \eqn{n_1, n_2, \theta_1, \theta_2}{n1, n2, theta1, theta2}. 
#'
#' @section Parameter recycling:
#'
#' It should be noted that the functions provided by sadists do \emph{not}
#' recycle their distribution parameters against the 
#' \code{x, p, q} or \code{n} parameters. This is in contrast to the 
#' common R idiom, and may cause some confusion. This is mostly for reasons
#' of performance, but also because some of the distributions have vector-valued
#' parameters; recycling over these would require the user to provide \emph{lists}
#' of parameters, which would be unpleasant.
#'
#' @section Legal Mumbo Jumbo:
#'
#' sadists is distributed in the hope that it will be useful,
#' but WITHOUT ANY WARRANTY; without even the implied warranty of
#' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#' GNU Lesser General Public License for more details.
#'
#' @template etc
#' @template ref-IP
#' @template ref-kprime
#' @template ref-Walck
#'
#' @import hypergeo PDQutils orthopolynom
#' @importFrom stats coefficients
#' @importFrom stats rchisq
#' @importFrom stats rnorm
#'
#' @note 
#' This package is maintained as a hobby. 
#'
#' @name sadists
#' @rdname sadists
#' @docType package
#' @title Some Additional Distributions
#' @keywords package
NULL

#' @title News for package 'sadists'
#'
#' @description
#' 
#' History of the 'sadists' package.
#'
#' \newcommand{\sadists}{\href{https://github.com/shabbychef/sadists}}
#'
#' @section \sadists{} Version 0.2.2 (2015-08-25) :
#' \itemize{
#' \item work around bad \code{rchisq} when \eqn{df=0=ncp}
#' }
#'
#' @section \sadists{} Version 0.2.1 (2015-06-12) :
#' \itemize{
#' \item shiny app (h/t Dean Attali).
#' }
#'
#' @section \sadists{} Version 0.2.0 (2015-04-01) :
#' \itemize{
#' \item add doubly non-central Beta and Eta distributions.
#' \item add (sum of) log chi-square distribution.
#' \item have products of chi-square depend on transform of this distribution.
#' }
#' @section \sadists{} Initial Version 0.1.0 (2015-03-07) :
#' \itemize{
#' \item first CRAN release.
#' }
#'
#' @name sadists-NEWS
#' @rdname NEWS
NULL


#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
