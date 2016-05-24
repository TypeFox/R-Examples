#' phenotypeData class
#'
#' Phenotype data
#'
#' @details Always a named list of two elements
#' \enumerate{
#' \item \code{y} numeric vector, phenotype
#' \item \code{yInfo} data.frame, additional information about observations
#' provied in .fam file
#' }
#'
#' @seealso \code{\link{read_phenotype}}
#' @name phenotypeData
NULL

#' Print phenotypeData class object
#'
#' @param x phenotypeData class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#'
#' @method print phenotypeData
print.phenotypeData <- function(x, ...){
  cat("Object of class phenotypeData\n")
  cat("$y: numeric vector of size", length(x$y), "\n")
  cat("$yInfo: data.frame\n")
  cat("\t", nrow(x$yInfo), " rows\n")
  cat("\t", ncol(x$yInfo), " columns\n")
}

#' Summary phenotypeData class object
#'
#' @param object phenotypeData class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#'
#' @method summary phenotypeData
summary.phenotypeData <- function(object, ...){
  cat("Object of class phenotypeData\n")
  cat("$y: numeric vector with phenotype of length", length(object$y), "\n")
  cat("$yInfo: matrix with additional information\n")
  cat("\t", nrow(object$yInfo), " observations\n")
  cat("\t", ncol(object$yInfo), " variables\n")
}
