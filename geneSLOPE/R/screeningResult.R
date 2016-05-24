#' screeningResult class
#'
#' A result of procedure for snp clumping produced by \code{\link{screen_snps}}
#'
#' @details Always a named list of eight elements
#' \enumerate{
#' \item \code{X} numeric matrix, consists of snps that passed screening
#' \item \code{y} numeric vector, phenotype
#' \item \code{X_info} data.frame, SNP info from .map file
#' \item \code{pVals} numeric vector, p-values from marginal tests for each snp
#' \item \code{numberOfSnps} numeric, total number of SNPs in .raw file
#' \item \code{selectedSnpsNumbers} numeric vector, which rows of \code{X_info}
#' data.frame are related to snps that passed screening
#' \item \code{pValMax} numeric, p-value used in screening procedure
#' \item \code{phenotypeInfo} data.frame, additional information about observations
#' provied in \code{\link{phenotypeData}} object
#' }
#'
#' @seealso \code{\link{phenotypeData}} \code{\link{screen_snps}}
#' @name screeningResult
NULL

#' Print function for class screeningResult class
#'
#' @param x screeningResult class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#'
#' @method print screeningResult
print.screeningResult <- function(x, ...){
  cat("Object of class screeningResult\n")
  cat("$X: numeric matrix\n")
  cat("\t", nrow(x$X), " rows\n")
  cat("\t", ncol(x$X), " columns\n")
  cat("$y: numeric phenotype vector of length", length(x$y), "\n")
  cat("$X_info: data.frame\n")
  cat("\t", nrow(x$X_info), " rows\n")
  cat("\t", ncol(x$X_info), " columns\n")
  cat("$pVals: numeric vector of length ", length(x$pVals), "\n")
  cat("$numberOfSnps: ", x$numberOfSnps, "\n")
  cat("$selectedSnpsNumbers: numeric vector of length" , length(x$selectedSnpsNumbers), "\n")
  cat("$pValMax: ", x$pValMax, "\n")
  cat("$phenotypeInfo: data.frame\n")
}


#' Summary function for class screeningResult
#'
#' @param object screeningResult class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#'
#' @method summary screeningResult
summary.screeningResult <- function(object, ...){
  cat("Object of class screeningResult\n")
  cat("$X: data matrix\n")
  cat("\t", nrow(object$X), " observations\n")
  cat("\t", ncol(object$X), " snps\n")
  cat(object$numberOfSnps, " SNPs were screened\n")
  cat(ncol(object$X), " snps had p-value smaller than ", object$pValMax,
      " in marginal test\n")
}
