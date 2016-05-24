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
#' Calculate the Distance Between Two Fuzzy Numbers
#'
#' @description
#' Currently, only Euclidean distance may be calculated.
#' We have \eqn{d_E^2(A,B) := \int_0^1  (A_L(\alpha)-B_L(\alpha))^2\,d\alpha,\int_0^1 + (A_U(\alpha)-B_U(\alpha))^2\,d\alpha
#' }{d_E^2(A,B) := int_0^1  (A_L(\alpha)-B_L(\alpha))^2 d\alpha + int_0^1  (A_U(\alpha)-B_U(\alpha))^2 d\alpha},
#' see (Grzegorzewski, 1988).
#'
#' @details
#' The calculation are done using numerical integration,
#'
#' @param e1 a fuzzy number
#' @param e2 a fuzzy number
#' @param type one of \code{"Euclidean"}, \code{"EuclideanSquared"}
#' @param ... additional arguments passed to \code{\link{integrate}}
#'
#' @return Returns the calculated distance, i.e. a single numeric value.
#'
#'
#' @exportMethod distance
#' @docType methods
#' @name distance
#' @family FuzzyNumber-method
#' @family DiscontinuousFuzzyNumber-method
#' @rdname distance-methods
#' @aliases distance,FuzzyNumber,FuzzyNumber-method
#'          distance,DiscontinuousFuzzyNumber,FuzzyNumber-method
#'          distance,FuzzyNumber,DiscontinuousFuzzyNumber-method
#'          distance,DiscontinuousFuzzyNumber,DiscontinuousFuzzyNumber-method
#'
#' @usage
#' \S4method{distance}{FuzzyNumber,FuzzyNumber}(e1, e2, type=c("Euclidean", "EuclideanSquared"), ...)
#'
#' \S4method{distance}{FuzzyNumber,DiscontinuousFuzzyNumber}(e1, e2, type=c("Euclidean", "EuclideanSquared"), ...)
#'
#' \S4method{distance}{DiscontinuousFuzzyNumber,FuzzyNumber}(e1, e2, type=c("Euclidean", "EuclideanSquared"), ...)
#'
#' \S4method{distance}{DiscontinuousFuzzyNumber,DiscontinuousFuzzyNumber}(e1, e2, type=c("Euclidean", "EuclideanSquared"), ...)
#'
#' @references
#' Grzegorzewski P., Metrics and orders in space of fuzzy numbers,
#' Fuzzy Sets and Systems 97, 1998, pp. 83-94.
setGeneric("distance",
           function(e1, e2, ...) standardGeneric("distance"))



setMethod(
   f="distance",
   signature(e1="FuzzyNumber", e2="FuzzyNumber"),
   definition=function(e1, e2,
                       type=c("Euclidean", "EuclideanSquared"), ...)
   {
      if (is.na(e1@lower(0)) || is.na(e2@lower(0))) return(NA_real_)
      type = match.arg(type)

      if (type == "Euclidean" || type == "EuclideanSquared")
      {
         dL <- integrate(function(alpha) {
            (  e1@a1+(e1@a2-e1@a1)*e1@lower(alpha)
               -e2@a1-(e2@a2-e2@a1)*e2@lower(alpha)
            )^2
         }, 0, 1, ...)$value

         dU <- integrate(function(alpha) {
            (  e1@a3+(e1@a4-e1@a3)*e1@upper(alpha)
               -e2@a3-(e2@a4-e2@a3)*e2@upper(alpha)
            )^2
         }, 0, 1, ...)$value

         if (type == "Euclidean") return (sqrt(dL+dU)) else return (dL+dU)
      } else
      {
         return(NA_real_)
      }
   }
)


setMethod(
   f="distance",
   signature(e1="DiscontinuousFuzzyNumber", e2="DiscontinuousFuzzyNumber"),
   definition=function(e1, e2,
      type=c("Euclidean", "EuclideanSquared"), ...)
   {
      if (is.na(e1@lower(0)) || is.na(e2@lower(0))) return(NA_real_)
      type = match.arg(type)

      if (type == "Euclidean" || type == "EuclideanSquared")
      {
         discontL <- c(e1@discontinuities.lower, e2@discontinuities.lower)
         discontL <- unique(sort(discontL))
         dL <- integrate_discont_val(function(alpha) {
            (  e1@a1+(e1@a2-e1@a1)*e1@lower(alpha)
              -e2@a1-(e2@a2-e2@a1)*e2@lower(alpha)
            )^2
         }, 0, 1, discontinuities=discontL, ...)

         discontU <- c(e1@discontinuities.upper, e2@discontinuities.upper)
         discontU <- unique(sort(discontU))
         dU <- integrate_discont_val(function(alpha) {
            (   e1@a3+(e1@a4-e1@a3)*e1@upper(alpha)
              - e2@a3-(e2@a4-e2@a3)*e2@upper(alpha)
            )^2
         }, 0, 1, discontinuities=discontU, ...)

         if (type == "Euclidean") return (sqrt(dL+dU)) else return (dL+dU)
      } else
      {
         return(NA_real_)
      }
   }
)




setMethod(
   f="distance",
   signature(e1="FuzzyNumber", e2="DiscontinuousFuzzyNumber"),
   definition=function(e1, e2,
      type=c("Euclidean", "EuclideanSquared"), ...)
   {
      return(distance(e2, e1, type=type, ...))
   }
)



setMethod(
   f="distance",
   signature(e1="DiscontinuousFuzzyNumber", e2="FuzzyNumber"),
   definition=function(e1, e2,
      type=c("Euclidean", "EuclideanSquared"), ...)
   {
      if (is.na(e1@lower(0)) || is.na(e2@lower(0))) return(NA_real_)
      type = match.arg(type)

      if (type == "Euclidean" || type == "EuclideanSquared")
      {
         dL <- integrate_discont_val(function(alpha) {
            (  e1@a1+(e1@a2-e1@a1)*e1@lower(alpha)
              -e2@a1-(e2@a2-e2@a1)*e2@lower(alpha)
            )^2
         }, 0, 1, discontinuities=e1@discontinuities.lower, ...)

         dU <- integrate_discont_val(function(alpha) {
            (  e1@a3+(e1@a4-e1@a3)*e1@upper(alpha)
              -e2@a3-(e2@a4-e2@a3)*e2@upper(alpha)
            )^2
         }, 0, 1, discontinuities=e1@discontinuities.upper, ...)

         if (type == "Euclidean") return (sqrt(dL+dU)) else return (dL+dU)
      } else
      {
         return(NA_real_)
      }
   }
)
