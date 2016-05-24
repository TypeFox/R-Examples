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
#' Print Basic Information on a Fuzzy Number
#'
#' @description
#' See \code{\link{as.character}} for more details.
#'
#' @details
#' The method \code{\link{as.character}} is called on given fuzzy number
#' object with default arguments.
#' The results are printed on \code{stdout}.
#'
#' @param object a fuzzy number
#' 
#' @return Does not return anything interesting.
#'
#' @usage
#' \S4method{show}{FuzzyNumber}(object)
#'
#' @exportMethod show
#' @name show
#' @aliases show,FuzzyNumber-method
#' @rdname show-methods
#' @family FuzzyNumber-method
#' @docType methods
if (!isGeneric("show"))
   setGeneric("show", function(object) standardGeneric("object"))



setMethod(
   f="show",
   signature(object="FuzzyNumber"),
   definition=function(object)
   {
      cat(as.character(object))
   }
)
