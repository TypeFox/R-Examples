# Copyright 2014-2014 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav
#
# This file is part of MarkowitzR.
#
# MarkowitzR is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MarkowitzR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MarkowitzR.  If not, see <http://www.gnu.org/licenses/>.

# Created: 2014.01.31
# Copyright: Steven E. Pav, 2014
# Author: Steven E. Pav
# Comments: Steven E. Pav

#' Inference on the Markowitz portfolio.
#' 
#' @section Markowitz Portfolio:
#'
#' Suppose \eqn{x} is a \eqn{p}-vector of returns of some assets with expected
#' value \eqn{\mu}{mu} and covariance \eqn{\Sigma}{Sigma}. The 
#' \emph{Markowitz Portfolio} is the portfolio 
#' \eqn{w = \Sigma^{-1}\mu}{w = Sigma^-1 mu}. Scale multiples of this portfolio
#' solve various portfolio optimization problems, among them
#' \deqn{\mathrm{argmax}_{w: w^{\top}\Sigma w \le R^2} \frac{\mu^{\top} w -
#' r_0}{\sqrt{w^{\top}\Sigma w}}}{argmax{ (mu'w - r0) / sqrt(w'Sigma w) :
#' w'Sigma w <= R^2}}
#'
#' This packages supports various statistical tests around the elements of 
#' the Markowitz Portfolio, and its Sharpe ratio, including the possibility of
#' hedging, and scalar conditional heteroskedasticity and conditional
#' expectation.
#'
#' @section Legal Mumbo Jumbo:
#'
#' MarkowitzR is distributed in the hope that it will be useful,
#' but WITHOUT ANY WARRANTY; without even the implied warranty of
#' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#' GNU Lesser General Public License for more details.
#'
#' @template etc
#' @template ref-SEP13
#' @template ref-MBJ99
#'
#' @references
#'
#' Bodnar, Taras and Okhrin, Yarema. "On the Product of Inverse Wishart
#' and Normal Distributions with Applications to Discriminant Analysis 
#' and Portfolio Theory." Scandinavian Journal of Statistics 38, no. 2 (2011):
#' 311--331. \url{http://dx.doi.org/10.1111/j.1467-9469.2011.00729.x}
#'
#' Markowitz, Harry. "Portfolio Selection." The Journal of Finance 7, no. 1
#' (1952): 77--91. \url{http://www.jstor.org/stable/2975974}
#'
#' Brandt, Michael W. "Portfolio Choice Problems." Handbook of Financial
#' Econometrics 1 (2009): 269--336. 
#' \url{https://faculty.fuqua.duke.edu/~mbrandt/papers/published/portreview.pdf}
#'
#' @import matrixcalc sandwich gtools
#'
#' @name MarkowitzR
#' @rdname MarkowitzR
#' @docType package
#' @title statistics concerning the Markowitz portfolio
#' @keywords package
#' 
#' This package is maintained as a hobby. 
#'
NULL

#' @title News for package 'MarkowitzR':
#'
#' \newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
#' \newcommand{\MarkowitzR}{\CRANpkg{MarkowitzR}}
#'
#' @section Changes in \MarkowitzR{} Version 0.1502 (2015-01-26) :
#' \itemize{
#' \item conform to CRAN rules.
#' }
#'
#' @section Changes in \MarkowitzR{} Version 0.1403 (2014-06-01) :
#' \itemize{
#' \item fix bug preventing multi-row hedging or constraint matrices.
#' }
#'
#' @section \MarkowitzR{} Initial Version 0.1402 (2014-02-14) :
#' \itemize{
#' \item first CRAN release.
#' }
#'
#' @name MarkowitzR-NEWS
#' @rdname NEWS
NULL

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
