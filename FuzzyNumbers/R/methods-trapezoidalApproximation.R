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
#' Trapezoidal Approximation of a Fuzzy Number
#'
#' @description
#' This method finds a trapezoidal approximation \eqn{T(A)}
#' of a given fuzzy number \eqn{A} by using the algorithm specified by the
#' \code{method} parameter.
#'
#' @details
#' \code{method} may be one of:
#' \enumerate{
#' \item \code{NearestEuclidean}: see (Ban, 2009); 
#' uses numerical integration, see \code{\link{integrateAlpha}}
#'
#' \item \code{Naive}:
#' We have core(A)==core(T(A)) and supp(A)==supp(T(A))
#'
#' \item \code{ExpectedIntervalPreserving}:
#' L2-nearest trapezoidal approximation preserving the expected interval given in
#' (Grzegorzewski, 2010; Ban, 2008; Yeh, 2008)
#' Unfortunately, for highly skewed membership functions 
#' this approximation operator may have
#' quite unfavourable behavior.
#' For example, if Val(A) < EV_{1/3}(A) or Val(A) > EV_{2/3}(A),
#' then it may happen that the core of the output
#' and the core of the original fuzzy number A are disjoint
#' (cf. Grzegorzewski, Pasternak-Winiarska, 2011)
#'
#' \item \code{SupportCoreRestricted}:
#' This method was proposed in (Grzegorzewski, Pasternak-Winiarska, 2011).
#' L2-nearest trapezoidal approximation with constraints
#' core(A) \eqn{\subseteq}{SUBSETS} core(T(A)) 
#' and supp(T(A)) \eqn{\subseteq}{SUBSETS} supp(A), i.e.
#' for which each point that surely belongs to A also belongs to T(A),
#' and each point that surely does not belong to A also does not belong to T(A).
#' }
#'
#' @usage
#' \S4method{trapezoidalApproximation}{FuzzyNumber}(object,
#'    method=c("NearestEuclidean", "ExpectedIntervalPreserving",
#'             "SupportCoreRestricted", "Naive"),
#'    ..., verbose=FALSE)
#'
#' @param object a fuzzy number
#' @param method character; one of: \code{"NearestEuclidean"} (default),
#' \code{"ExpectedIntervalPreserving"},
#' \code{"SupportCoreRestricted"},
#' \code{"Naive"}
#' @param verbose logical; should some technical details on the computations being performed be printed?
#' @param ... further arguments passed to \code{\link{integrateAlpha}}
#' 
#' @return Returns a \code{\link{TrapezoidalFuzzyNumber}} object.
#'
#' @exportMethod trapezoidalApproximation
#' @name trapezoidalApproximation
#' @aliases trapezoidalApproximation,FuzzyNumber-method
#' @rdname trapezoidalApproximation-methods
#' @docType methods
#' @family approximation
#' @family FuzzyNumber-method
#' @references
#' Ban A.I. (2008), Approximation of fuzzy numbers by trapezoidal fuzzy numbers
#' preserving the expected interval, Fuzzy Sets and Systems 159, pp. 1327-1344.
#' 
#' Ban A.I. (2009), On the nearest parametric approximation of a fuzzy number 
#' - Revisited, Fuzzy Sets and Systems 160, pp. 3027-3047.
#' 
#' Grzegorzewski P. (2010), Algorithms for trapezoidal 
#' approximations of fuzzy numbers
#' preserving the expected interval, In: Bouchon-Meunier B. et al (Eds.),
#' Foundations of Reasoning Under Uncertainty, Springer, pp. 85-98.
#' 
#' Grzegorzewski P, Pasternak-Winiarska K. (2011), Trapezoidal
#'  approximations of fuzzy numbers
#' with restrictions on the support and core, Proc. EUSFLAT/LFA 2011,
#'  Atlantic Press, pp. 749-756.
#' 
#' Yeh C.-T. (2008), Trapezoidal and triangular approximations 
#' preserving the expected interval,
#' Fuzzy Sets and Systems 159, pp. 1345-1353.
#'
#' @examples
#' (A <- FuzzyNumber(-1, 0, 1, 40,
#'    lower=function(x) sqrt(x), upper=function(x) 1-sqrt(x)))
#' (TA <- trapezoidalApproximation(A,
#'    "ExpectedIntervalPreserving")) # Note that the cores are disjoint!
#' expectedInterval(A)
#' expectedInterval(TA)
setGeneric("trapezoidalApproximation",
           function(object, ...)
              standardGeneric("trapezoidalApproximation"));


setMethod(
   f="trapezoidalApproximation",
   signature(object="FuzzyNumber"),
   definition=function(
      object,
      method=c("NearestEuclidean", "ExpectedIntervalPreserving", "SupportCoreRestricted", "Naive"),
      ...,
      verbose=FALSE)
   {
      method <- match.arg(method)

      if (method == "Naive")
      {
         return(trapezoidalApproximation_Naive(object))
      }
      else if (method == "NearestEuclidean")
      {
         return(trapezoidalApproximation_NearestEuclidean(object, ..., verbose=verbose))
      }
      else if (method == "ExpectedIntervalPreserving")
      {
         return(trapezoidalApproximation_ExpectedIntervalPreserving(object, ..., verbose=verbose))
      }
      else if (method == "SupportCoreRestricted")
      {
         return(trapezoidalApproximation_SupportCoreRestricted(object, ..., verbose=verbose))
      }
   }
)


# internal function
trapezoidalApproximation_Naive <- function(object)
{
   return(TrapezoidalFuzzyNumber(object@a1, object@a2, object@a3, object@a4))
}


# internal function
trapezoidalApproximation_NearestEuclidean <- function(object, ..., verbose)
{
   # calculate integrals:
   if (is.na(object@lower(0)) || is.na(object@upper(0)))
      stop("Integral of alphacut bounds cannot be computed")

   expected.interval <- expectedInterval(object, ...)
   alpha.interval <- alphaInterval(object, ...)

   intLower <- expected.interval[1]
   intUpper <- expected.interval[2]
   intAlphaTimesLower <- alpha.interval[1]
   intAlphaTimesUpper <- alpha.interval[2]

   # Here we use the method given in (Ban, 2009)

   if (intAlphaTimesUpper-intAlphaTimesLower >= (intUpper-intLower)/3 )
      # i.e. if ambiguity(A) >= width(A)/3
   { # (i)
      if (verbose) cat("Using Case (i) of (Corollary 8; Ban, 2009)\n")
      a1    <-  4*intLower-6*intAlphaTimesLower
      a2    <- -2*intLower+6*intAlphaTimesLower
      a3    <- -2*intUpper+6*intAlphaTimesUpper
      a4    <-  4*intUpper-6*intAlphaTimesUpper
   } else
      if (-intLower+3*intAlphaTimesLower-3*intUpper+5*intAlphaTimesUpper > 0)
      { # (iii)
         if (verbose) cat("Using Case (iii) of (Corollary 8; Ban, 2009)\n")

         a1 <-             (16*intLower-18*intAlphaTimesLower-2*intUpper)/5
         a2 <- a3 <- a4 <- (-2*intLower+ 6*intAlphaTimesLower+4*intUpper)/5
      } else
         if (3*intLower-5*intAlphaTimesLower+intUpper-3*intAlphaTimesUpper > 0)
         { # (iv)
            if (verbose) cat("Using Case (iv) of (Corollary 8; Ban, 2009)\n")
            a1 <- a2 <- a3 <- ( 4*intLower- 2*intUpper+ 6*intAlphaTimesUpper)/5
            a4 <-             (-2*intLower+16*intUpper-18*intAlphaTimesUpper)/5
         } else
         { # (ii)
            if (verbose) cat("Using Case (ii) of (Corollary 8; Ban, 2009)\n")

            a1       <- ( 7*intLower-9*intAlphaTimesLower+1*intUpper-3*intAlphaTimesUpper)/2
            a2 <- a3 <- (-2*intLower+6*intAlphaTimesLower-2*intUpper+6*intAlphaTimesUpper)/2
            a4       <- ( 1*intLower-3*intAlphaTimesLower+7*intUpper-9*intAlphaTimesUpper)/2
         }

   return(TrapezoidalFuzzyNumber(a1,a2,a3,a4))
}


# internal function
trapezoidalApproximation_ExpectedIntervalPreserving <- function(object, ..., verbose)
{
   # calculate integrals:
   if (is.na(object@lower(0)) || is.na(object@upper(0)))
      stop("Integral of alphacut bounds cannot be computed")

   expected.interval <- expectedInterval(object, ...)
   alpha.interval <- alphaInterval(object, ...)

   intLower <- expected.interval[1]
   intUpper <- expected.interval[2]
   intAlphaTimesLower <- alpha.interval[1]
   intAlphaTimesUpper <- alpha.interval[2]

   # Here we use the method given in (Grzegorzewski, 2010)

   if (intAlphaTimesUpper-intAlphaTimesLower >= (intUpper-intLower)/3 )
      # i.e. if ambiguity(A) >= width(A)/3
   {
      a1 <-  4*intLower-6*intAlphaTimesLower
      a2 <- -2*intLower+6*intAlphaTimesLower
      a3 <- -2*intUpper+6*intAlphaTimesUpper
      a4 <-  4*intUpper-6*intAlphaTimesUpper
   } else
   {
      Eval13 <- 2*intLower/3 +   intUpper/3 # Weighted Expected Value w=1/3
      Eval23 <-   intLower/3 + 2*intUpper/3 # Weighted Expected Value w=2/3
      Val <- intAlphaTimesLower + intAlphaTimesUpper # Value

      if (Eval13 <= Val && Val <= Eval23)
      {
         a1 <-       3*intLower +   intUpper - 3*intAlphaTimesLower - 3*intAlphaTimesUpper
         a2 <- a3 <-  -intLower -   intUpper + 3*intAlphaTimesLower + 3*intAlphaTimesUpper
         a4 <-         intLower + 3*intUpper - 3*intAlphaTimesLower - 3*intAlphaTimesUpper
      } else if (Val < Eval13)
      {
         a1 <- a2 <- a3 <- intLower
         a4 <- 2*intUpper - intLower
      } else
      {
         a1 <- 2*intLower - intUpper
         a2 <- a3 <- a4 <- intUpper
      }
   }

   return(TrapezoidalFuzzyNumber(a1, a2, a3, a4))

}



# internal function
trapezoidalApproximation_SupportCoreRestricted <- function(object, ..., verbose)
{
   # calculate integrals:
   if (is.na(object@lower(0)) || is.na(object@upper(0)))
      stop("Integral of alphacut bounds cannot be computed")

   expected.interval <- expectedInterval(object, ...)
   alpha.interval <- alphaInterval(object, ...)

   intLower <- expected.interval[1]
   intUpper <- expected.interval[2]
   intAlphaTimesLower <- alpha.interval[1]
   intAlphaTimesUpper <- alpha.interval[2]

   # Here we use the method given in (Grzegorzewski, Pasternak-Winiarska, 2011)

   u1 <- 4*intLower - 6*intAlphaTimesLower
   u2 <- 6*intAlphaTimesLower - 2*intLower
   if (object@a2 > u2)
   {
      if (object@a1 < u1)
      {
         a1 <- u1
         a2 <- u2
      } else
      {
         a1 <- object@a1
         a2 <- 2*intLower-object@a1
      }
   } else
   {
      if (object@a1 < u1)
      {
         a1 <- 2*intLower-object@a2
         a2 <- object@a2
      } else
      {
         a1 <- object@a1
         a2 <- object@a2
      }
   }


   u3 <- 6*intAlphaTimesUpper - 2*intUpper
   u4 <- 4*intUpper - 6*intAlphaTimesUpper
   if (object@a4 > u4)
   {
      if (object@a3 < u3)
      {
         a3 <- u3
         a4 <- u4
      } else
      {
         a3 <- object@a3
         a4 <- 2*intUpper-object@a3
      }
   } else
   {
      if (object@a3 < u3)
      {
         a3 <- 2*intUpper-object@a4
         a4 <- object@a4
      } else
      {
         a3 <- object@a3
         a4 <- object@a4
      }
   }

   return(TrapezoidalFuzzyNumber(a1, a2, a3, a4))
}
