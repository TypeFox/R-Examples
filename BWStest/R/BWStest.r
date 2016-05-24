# Copyright 2016-2016 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav

# This file is part of BWStest.
#
# BWStest is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BWStest is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with BWStest.  If not, see <http://www.gnu.org/licenses/>.


#' Baumgartner Weiss Schindler test.
#' 
#' @section Background:
#'
#' The Baumgartner Weiss Schindler test is a two sample test of the null 
#' that the samples come from the same probability distribution, similar
#' to the Kolmogorv-Smirnov, Wilcoxon, and Cramer-Von Mises tests. It is
#' similar to the Cramer-Von Mises test in that it estimates the
#' square norm of the difference in CDFs of the two samples. However, the
#' Baumgartner Weiss Schindler test weights the integral by the variance
#' of the difference in CDFs, "[emphasizing] the tails of the distributions,
#' which increases the power of the test for a lot of applications."
#'
#' @section Legal Mumbo Jumbo:
#'
#' BWStest is distributed in the hope that it will be useful,
#' but WITHOUT ANY WARRANTY; without even the implied warranty of
#' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#' GNU Lesser General Public License for more details.
#'
#' @template etc
#' @template ref-bws
#' @template ref-modtests
#' @name BWStest-package
#' @rdname BWStest-package
#' @docType package
#' @title Baumgartner Weiss Schindler test of equal distributions.
#' @keywords package
#' @import Rcpp 
#' @useDynLib BWStest
#' @importFrom Rcpp evalCpp
#' @importFrom stats ecdf
#' @importFrom memoise memoise
#' @exportPattern "^[[:alpha:]]+"
#'
NULL

#' @title News for package 'BWStest':
#'
#' @description 
#'
#' News for package 'BWStest'
#'
#' \newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
#' \newcommand{\cranBWStest}{\CRANpkg{BWStest}}
#' \newcommand{\BWStest}{\href{https://github.com/shabbychef/BWStest}}
#'
#' @section \BWStest{} Initial Version 0.2.0 (2016-04-29) :
#' \itemize{
#' \item Adding Murakami statistics.
#' }
#'
#' @section \BWStest{} Initial Version 0.1.0 (2016-04-07) :
#' \itemize{
#' \item First CRAN release.
#' }
#'
#' @section \BWStest{} Initial Version 0.0.0 (2016-04-06) :
#' \itemize{
#' \item Start work
#' }
#'
#' @name BWStest-NEWS
#' @rdname NEWS
NULL

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
