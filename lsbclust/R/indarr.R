#' Create Array of Indicator Matrices
#' 
#' This function takes a \code{matrix} or \code{data.frame} and the number of rating categories \code{maxcat}
#' and produces a three-way array of \code{m} by \code{maxcat} indicator matrices, one for each of the \code{n} rows.
#' The input \code{x} must be a \code{matrix} or \code{data.frame} of dimensions \code{n} by \code{m} 
#' which contains the ratings on a scale of 1 to \code{maxcat} for \code{m} items. Note that missing values
#' (\code{NA}'s) will not appear in the columns.
#' 
#' @param x a \code{matrix} of \code{data.frame}
#' @param maxcat an integer indicating the maximum of the rating scale (which is assumed to start with 1)
#' @param na.add logical indicating whether to add a designated category for missings or not. Defaults to TRUE.
#' @author Pieter C. Schoonees
#' @export
#' @return A list of rating by item indicator matrices.
#' @examples
#' data("lov")
#' arr <- indarr(lov[1:10, 1:9], maxcat = 9)
#' str(arr)
indarr <- function(x, maxcat, na.add = TRUE){
  stopifnot(class(x)[1] %in% c("matrix", "data.frame"))
  if(any(!complete.cases(x)) && na.add) {
    rnames <- c(1:maxcat, "NA")
    maxcat <- maxcat + 1
    x[is.na(x)] <- maxcat
  } else rnames <- 1:maxcat
  cnames <- colnames(x)
  x <- as.matrix(x)
  dm <- diag(maxcat)
  afun <- function(y, dm){
    tmp <- dm[, y]
    return(tmp)
  }
  out <- array(apply(x, MARGIN = 1, FUN =  afun, dm = dm), dim = c(maxcat, ncol(x), nrow(x)))
  out[is.na(out)] <- 0
  rownames(out) <- rnames
  colnames(out) <- cnames
  return(out)
}