# Copyright 2015-2015 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav
#
# This file is part of madness.
#
# madness is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# madness is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with madness.  If not, see <http://www.gnu.org/licenses/>.

# Created: 2015.11.16
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav
# Comments: Steven E. Pav

#' Automatic Differentiation of Matrix Operations.
#'
#' @section Legal Mumbo Jumbo:
#'
#' madness is distributed in the hope that it will be useful,
#' but WITHOUT ANY WARRANTY; without even the implied warranty of
#' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#' GNU Lesser General Public License for more details.
#'
#' @template etc
#'
#' @references
#'
#' Griewank, Andreas and Walther, Andrea. "Evaluating Derivatives: principles and techniques of algorithmic differentiation."
#' SIAM (2008).
#'
#' Petersen, Kaare Brandt and Pedersen, Michael Syskind. "The Matrix Cookbook."
#' Technical University of Denmark (2012). 
#' \url{http://www2.imm.dtu.dk/pubdb/p.php?3274}
#'
#' Magnus, Jan R. and Neudecker, H. "Matrix Differential Calculus with Applications in Statistics and Econometrics."
#' 3rd Edition. Wiley Series in Probability and Statistics: Texts and References Section (2007).
#' \url{http://www.janmagnus.nl/misc/mdc2007-3rdedition}
#'
#' Magnus, Jan R. and Neudecker, H. "The elimination matrix: some lemmas and applications," 
#' SIAM Journal on Algebraic Discrete Methods 1, no. 4 (1980): 422-449.
#' \url{http://www.janmagnus.nl/papers/JRM008.pdf}
#'
#' Magnus, Jan R. and Neudecker, H. "Symmetry, 0-1 Matrices and Jacobians,"
#' Econometric Theory 2 (1986): 157-190.
#' \url{http://www.janmagnus.nl/papers/JRM014.pdf},
#'
#' Fackler, Paul L. "Notes on Matrix Calculus." (2005).
#' \url{http://www4.ncsu.edu/~pfackler/MatCalc.pdf}
#'
#' @import matrixcalc expm methods
#' @importFrom methods cbind2
#' @importFrom methods rbind2
#' @importFrom methods coerce
#' @importFrom methods show
#' @importFrom methods kronecker
#' @importFrom stats coef
#' @importFrom stats lm
#' @importFrom stats na.omit
#' @importFrom stats vcov
#' @importFrom utils head
#'
#' @name madness-pkg
#' @rdname madness-pkg
#' @docType package
#' @title Multivariate Automatic Differentiation.
#' @keywords package
#' @template etc
#' @note
#' 
#' This package is maintained as a hobby. 
#'
NULL

#' @title News for package \sQuote{madness}:
#'
#' @description
#'
#' News for package \sQuote{madness}.
#'
#' \newcommand{\pkg}{#1}
#' \newcommand{\CRANpkg}{\href{https://CRAN.R-project.org/package=#1}{\pkg{#1}}}
#' \newcommand{\madness}{\CRANpkg{madness}}
#'
#' @section \madness{} Version 0.1.0.400 (2016-01-12) :
#' \itemize{
#' \item adding \code{max} and \code{min}.
#' }
#'
#' @section \madness{} Version 0.1.0.300 (2016-01-10) :
#' \itemize{
#' \item adding \code{eigen}.
#' }
#'
#' @section \madness{} Version 0.1.0.200 (2016-01-07) :
#' \itemize{
#' \item exporting \code{diag}.
#' }
#'
#' @section \madness{} Version 0.1.0 (2015-12-15) :
#' \itemize{
#' \item first CRAN release.
#' }
#'
#' @section \madness{} Initial Version 0.0.0.5000 (2015-12-01) :
#' \itemize{
#' \item first github release.
#' }
#'
#' @name madness-NEWS
#' @rdname NEWS
NULL


# 2FIX
#
# expm logm (these will require Matrix package ...)
# eigs?
#
# %o%   (note this is a *huge* outer product potentially)
#

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
