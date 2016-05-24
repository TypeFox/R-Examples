#' Summarizing ALK Estimations
#' 
#' \code{summary} method for class "\code{ALKr}". Returns the mean
#' length-at-age, the variance of the length-at-age, and the proportion of the
#' age class in the population, for each age class.
#' 
#' @param object An \code{ALKr} object.
#' @param length_classes A vector with the length value to be used for each
#' length class. Defaults to the row names of \code{x}, which must be coercible
#' as numerics.
#' @param ... other arguments (currently ignored)
#' 
#' @return A \code{summary_ALKr} object, containing the mean length-at-age, the
#' variance of the length-at-age, and the proportion of the age class in the
#' population, for each age class. The name of the method used to calculate the
#' ALK and its parameters are also included.
#' 
#' @examples
#' data(hom)
#' cALK <- classic_ALK(hom$otoliths[[1]], fi = hom$F1992)
#' summary(cALK)
#' hhALK <- hoenig_heisey(hom$otoliths[[1]], fi1 = hom$F1992, fi2 = hom$F1993)
#' summary(hhALK)
#' 
#' @method summary ALKr
#' @export
summary.ALKr <- function(object, length_classes = NULL, ...) {
  
  if (is.null(length_classes)) {
    length_classes <- as.numeric(rownames(object@alk))  
  }
  
  nj <- colSums(object@N)
  i <- colSums(length_classes * object@N)
  lj <- i / nj
  
  vlj <- colSums(object@N * length_classes^2) / (nj - 1) - 2 * lj * i / (nj - 1) + lj^2 * nj / (nj - 1)
  vlj[vlj <= 0] <- NaN
  
  new("summary_ALKr",
      pj = nj / sum(nj),
      mean_lj = lj,
      var_lj = vlj,
      method = object@method,
      parameters = object@parameters,
      name = object@name,
      description = object@description
  )
}
