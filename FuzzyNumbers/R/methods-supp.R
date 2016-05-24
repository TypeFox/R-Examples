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
#' Calculate the Support of a Fuzzy Number
#'
#' @description
#' We have \eqn{\mathrm{supp}(A) := [a1,a4]}{supp(A) := [a1,a4]}.
#' This gives the values that a fuzzy number possibly may represent.
#'
#' @param object a fuzzy number
#' 
#' @return Returns a numeric vector of length 2.
#'
#' @exportMethod supp
#' @docType methods
#' @name supp
#' @family FuzzyNumber-method
#' @family alpha_cuts
#' @rdname supp-methods
#' @aliases supp,FuzzyNumber-method
#' @usage
#' \S4method{supp}{FuzzyNumber}(object)
setGeneric("supp",
           function(object) standardGeneric("supp"))



setMethod(
   f="supp",
   signature(object="FuzzyNumber"),
   definition=function(object)
   {
      c(object@a1, object@a4)
   }
)
