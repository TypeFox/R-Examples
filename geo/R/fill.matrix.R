#' Relace elements of a matrix
#' 
#' Replace (fill) elements of a matrix (or data.frame) with a value for given
#' pairs of row and column indices.
#' 
#' 
#' @param outcome Input matrix/data.frame
#' @param x Value or values to replace/fill with
#' @param rownr Row index/indices
#' @param dalknr Column index/indices
#' @return Matrix or data.frame with given values replaced.
#' @note Probably redundant, not called by any geo-function, the same effect
#' could be achieved with an assignment to a matrix with an index-matrix of the
#' values in rownr and dalknr: \code{mat[matrix(c(rownr, dalknr), ncol = 2)] <-
#' x}
#' @keywords manip
#' @export fill.matrix
fill.matrix <-
function(outcome, x, rownr, dalknr)
{
        ind <- nrow(outcome) * (dalknr - 1) + rownr
        outcome[ind] <- x
        return(outcome)
}

