## This file is part of the FuzzyNumbers library.
##
## Copyright 2012-2014 Marek Gagolewski
##
##
## FuzzyNumbers is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## FuzzyNumbers is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with FuzzyNumbers. If not, see <http://www.gnu.org/licenses/>.


#' @title Tools to Deal with Fuzzy Numbers
#'
#' @description
#' \pkg{FuzzyNumbers} is an open source (LGPL 3) package for R.
#' It provides S4 classes and methods to deal with Fuzzy Numbers.
#' The package may be used by the practitioners as well as by the researchers
#' in fuzzy numbers theory (e.g. for testing new algorithms,
#' generating numerical examples, preparing figures).
#'
#' @details
#' Fuzzy set theory lets us quite intuitively represent
#' imprecise or vague information. Fuzzy numbers, which form a particular
#' subclass of fuzzy sets of the real line, play a significant role
#' in many important theoretical and/or practical considerations.
#' This is because we often describe our knowledge about objects
#' through numbers, e.g. "I'm about 180 cm tall"
#' or "The rocket was launched between 2 and 3 p.m.".
#'
#' For the formal definition of a fuzzy number
#' please refer to the \code{\linkS4class{FuzzyNumber}} man page.
#' Note that this package also deals with particular types
#' of fuzzy numbers like trapezoidal, piecewise linear, or ``parametric'' FNs
#' (see \code{\linkS4class{TrapezoidalFuzzyNumber}}
#' \code{\linkS4class{PiecewiseLinearFuzzyNumber}},
#' \code{\linkS4class{PowerFuzzyNumber}}, and *EXPERIMENTAL*
#' \code{\linkS4class{DiscontinuousFuzzyNumber}})
#'
#' The package aims to provide the following functionality:
#' \enumerate{
#'    \item Representation of arbitrary fuzzy numbers
#' (including FNs with discontinuous side functions and/or alpha-cuts),
#' as well as their particular types, e.g. trapezoidal and piecewise linear fuzzy numbers,
#'    \item Defuzzification and approximation by triangular
#'      and piecewise linear FNs (see e.g. \code{\link{expectedValue}},
#'      \code{\link{value}}, \code{\link{trapezoidalApproximation}},
#'      \code{\link{piecewiseLinearApproximation}}),
#'    \item Visualization of FNs (see \code{\link{plot}}, \code{\link{as.character}}),
#'    \item Basic operations on FNs (see e.g. \code{\link{fapply}} and \link{Arithmetic}),
#'    \item Aggregation of FNs **TO DO**,
#'    \item Ranking of FNs **TO DO**,
#'    \item Random FN generation **TO DO**,
#'    \item \dots
#' }
#'
#' Please feel free to send any comments and feature requests to the author
#' (see his homepage at \url{http://gagolewski.rexamine.com/}).
#'
#' For a complete list of classes and methods
#' call \code{help(package="FuzzyNumbers")}.
#' Moreover, you will surely be interested in a step-by-step guide to the
#' package usage and features which is available at the project's webpage.
#' \cr\cr
#'
#'
#'
#' \bold{Keywords}: Fuzzy Numbers, Fuzzy Sets, Shadowed Sets,
#' Trapezoidal Approximation, Piecewise Linear Approximation,
#' Approximate Reasoning, Imprecision, Vagueness, Randomness.
#'
#' \bold{Acknowledgments}: Many thanks to Jan Caha,
#' Przemyslaw Grzegorzewski, Lucian Coroianu, and Pablo Villacorta Iglesias
#'  for stimulating discussion.
#'
#'  The development of the package in March-June 2013 was partially supported
#'  by the European Union from resources of the European Social Fund, Project PO KL
#'  ``Information technologies: Research and their interdisciplinary
#'  applications'', agreement UDA-POKL.04.01.01-00-051/10-00.
#'
#' @name FuzzyNumbers-package
#' @docType package
#' @import methods
#' @import grDevices
#' @import graphics
#' @import stats
#' @author Marek Gagolewski \email{gagolews@@rexamine.com},\cr
#'  with contributions from Jan Caha
#' @references
#' \pkg{FuzzyNumbers} Homepage, \url{http://FuzzyNumbers.rexamine.com/}.
#' 
#' Ban A.I. (2008), Approximation of fuzzy numbers by trapezoidal fuzzy numbers
#' preserving the expected interval, Fuzzy Sets and Systems 159, pp. 1327-1344.
#' 
#' Ban A.I. (2009), On the nearest parametric approximation of a fuzzy number - Revisited,
#' Fuzzy Sets and Systems 160, pp. 3027-3047.
#' 
#' Bodjanova S. (2005), Median value and median interval of a fuzzy number,
#' Information Sciences 172, pp. 73-89.
#' 
#' Chanas S. (2001), On the interval approximation of a fuzzy number,
#' Fuzzy Sets and Systems 122, pp. 353-356.
#' 
#' Coroianu L., Gagolewski M., Grzegorzewski P. (2013),
#' Nearest Piecewise Linear Approximation of Fuzzy Numbers,
#' Fuzzy Sets and Systems 233, pp. 26-51.
#' 
#' Coroianu L., Gagolewski M., Grzegorzewski P.,
#' Adabitabar Firozja M., Houlari T. (2014),
#' Piecewise linear approximation of fuzzy numbers preserving 
#' the support and core, In: Laurent A. et al. (Eds.), 
#' Information Processing and Management of Uncertainty in 
#' Knowledge-Based Systems, Part II (CCIS 443), Springer, pp. 244-254.
#' 
#' Delgado M., Vila M.A., Voxman W. (1998), 
#' On a canonical representation of a fuzzy number,
#' Fuzzy Sets and Systems 93, pp. 125-135.
#' 
#' Dubois D., Prade H. (1978), Operations on fuzzy numbers, 
#' Int. J. Syst. Sci. 9, pp. 613-626.
#' 
#' Dubois D., Prade H. (1987a), The mean value of a fuzzy number, 
#' Fuzzy Sets and Systems 24, pp. 279-300.
#' 
#' Dubois D., Prade H. (1987b), Fuzzy numbers: An overview, In: Analysis of Fuzzy
#' Information. Mathematical Logic, vol. I, CRC Press, pp. 3-39.
#' 
#' Grzegorzewski P. (2010), Algorithms for trapezoidal approximations of fuzzy numbers
#' preserving the expected interval, In: Bouchon-Meunier B. et al (Eds.),
#' Foundations of Reasoning Under Uncertainty, Springer, pp. 85-98.
#' 
#' Grzegorzewski P. (1998), Metrics and orders in space of fuzzy numbers,
#' Fuzzy Sets and Systems 97, pp. 83-94.
#' 
#' Grzegorzewski P,. Pasternak-Winiarska K. (2011), 
#' Trapezoidal approximations of fuzzy numbers
#' with restrictions on the support and core,
#'  Proc. EUSFLAT/LFA 2011, Atlantic Press, pp. 749-756.
#' 
#' Klir G.J., Yuan B. (1995), Fuzzy sets and fuzzy logic. 
#' Theory and applications, Prentice Hall, New Jersey.
#' 
#' Stefanini L., Sorini L. (2009), Fuzzy arithmetic with 
#' parametric LR fuzzy numbers,
#' In: Proc. IFSA/EUSFLAT 2009, pp. 600-605.
#' 
#' Yeh C.-T. (2008), Trapezoidal and triangular approximations 
#' preserving the expected interval,
#' Fuzzy Sets and Systems 159, pp. 1345-1353.
#' 
#' @import methods
#' @import grDevices
#' @import graphics
#' @import stats
invisible(NULL)
