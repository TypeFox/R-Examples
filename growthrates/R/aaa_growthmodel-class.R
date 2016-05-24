#' Class of Growth Model Functions
#'
#' This class is used for the parametric \code{grow_\dots} functions of the
#' package and can also be used for user-defined functions to describe
#' time-dependent growth of organisms.
#'
#' @seealso the constructor function \code{\link{growthmodel}} how to create
#'   instances of this class.
#' @name growthmodel-class
#' @exportClass growthmodel
#'
setOldClass("growthmodel") # S3 class

#' Union Class of Growth Model or Function
#'
#' This class union comprises parametric model functions from class
#' \code{growthmodel} and ordinary functions to describe time-dependent
#' growth of organisms.
#'
#' @seealso the constructor function \code{\link{growthmodel}} how to create
#'   instances of class \code{growthmodel}.
#'
#' @name function_growthmodel-class
#' @exportClass function_growthmodel
#'
setClassUnion("function_growthmodel", c("growthmodel", "function"))

