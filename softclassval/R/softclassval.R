##' Soft classification performance measures
##'
##' Extension of sensitivity, specificity, positive and negative predictive value to continuous
##' predicted and reference memberships in [0, 1].
##' @name softclassval-package
##' @docType package
##' @author C. Beleites
##' @import arrayhelpers
{
  if (!requireNamespace ("svUnit", quietly = TRUE)){
    `.test<-` <- function (f, value) {
      class (value) <-  c ("svTest", "function")
      attr (f, "test") <- value
      f
    }
  } else {
    `.test<-` <- svUnit::`test<-`
  }
  
  checkEqualsOrdered <- function (target, current, ...)
    checkEquals (target [order (names (target))], current [order (names (current))], ...)

  checkEqualAttributes <- function (target, current, ...)
    checkEqualsOrdered (attributes (target), attributes (current), ...) # TODO: exclusion list
}

