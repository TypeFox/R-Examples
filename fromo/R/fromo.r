# Copyright 2016-2016 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav

# This file is part of fromo.
#
# fromo is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# fromo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with fromo.  If not, see <http://www.gnu.org/licenses/>.


#' Fast Robust Moments.
#' 
#' @section Robust Moments:
#'
#' Welford described a method for 'robust' one-pass computation of the
#' standard deviation. By 'robust', we mean robust to round-off caused
#' by a large shift in the mean. This method was generalized by Terriberry,
#' and Bennett \emph{et. al.} to the case of higher-order moments. 
#' This package provides those algorithms for computing moments.
#'
#' Generally we should find that the stock implementations of \code{sd},
#' \code{skewness} and so on are \emph{already} robust and likely using
#' these algorithms under the hood. This package was written for a few
#' reasons:
#' \enumerate{
#' \item As an exercise to learn Rcpp.
#' \item Often I found I needed the first \eqn{k} moments. For example,
#' when computing the Z-score, the standard deviation and mean must be
#' computed separately, which is inefficient. Similarly Merten's correction
#' for the standard error of the Sharpe ratio uses the first four moments.
#' These are all computed as a side effect of computation of the kurtosis,
#' but discarded by the standard methods.
#' }
#'
#' @section Legal Mumbo Jumbo:
#'
#' fromo is distributed in the hope that it will be useful,
#' but WITHOUT ANY WARRANTY; without even the implied warranty of
#' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#' GNU Lesser General Public License for more details.
#'
#' @template etc
#' @template ref-romo
#' @note 
#'
#' This package was developed as an exercise in learning Rcpp.
#'
#' @name fromo-package
#' @rdname fromo-package
#' @docType package
#' @title fast robust moments
#' @keywords package
#' @import Rcpp methods
#' @useDynLib fromo
#' @importFrom Rcpp evalCpp
#' @importFrom methods new 
#' @exportPattern "^[[:alpha:]]+"
#'
NULL

#' @title News for package 'fromo':
#'
#' @description 
#'
#' News for package 'fromo'
#'
#' \newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
#' \newcommand{\cranfromo}{\CRANpkg{fromo}}
#' \newcommand{\fromo}{\href{https://github.com/shabbychef/fromo}}
#'
#' @section \fromo{} Version 0.1.3 (2016-04-04) :
#' \itemize{
#' \item submit to CRAN
#' }
#'
#' @section \fromo{} Initial Version 0.1.0 (2016-03-25) :
#' \itemize{
#' \item start work
#' }
#'
#' @name fromo-NEWS
#' @rdname NEWS
NULL

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
