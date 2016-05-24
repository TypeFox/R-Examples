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



#' @title
#' Calculate the Core of a Fuzzy Number
#'
#' @description
#' We have \eqn{\mathrm{core}(A) := [a2,a3]}{core(A) := [a2,a3]}.
#' This gives the values that a fuzzy number necessarily represents.
#'
#' @param object a fuzzy number
#' 
#' @return Returns a numeric vector of length 2.
#'
#' @exportMethod core
#' @docType methods
#' @name core
#' @family FuzzyNumber-method
#' @family alpha_cuts
#' @rdname core-methods
#' @aliases core,FuzzyNumber-method
#' @usage
#' \S4method{core}{FuzzyNumber}(object)
setGeneric("core",
           function(object) standardGeneric("core"))



setMethod(
   f="core",
   signature(object="FuzzyNumber"),
   definition=function(object)
   {
      c(object@a2, object@a3)
   }
)
