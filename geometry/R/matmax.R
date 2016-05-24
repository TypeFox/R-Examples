##' Row-wise matrix functions
##' 
##' Compute maximum or minimum of each row, or sort each row of a matrix, or a
##' set of (equal length) vectors.
##' 
##' 
##' @aliases matmax matmin matsort matorder
##' @param \dots A numeric matrix or a set of numeric vectors (that are
##' column-wise bind together into a matrix with \code{cbind}).
##' @return \code{matmin} and \code{matmax} return a vector of length
##' \code{nrow(cbind(...))}.  \code{matsort} returns a matrix of dimension
##' \code{dim(cbind(...))} with in each row of \code{cbind(...)} sorted.
##' \code{matsort(x)} is a lot faster than, e.g., `\code{t(apply(x,1,sort))}',
##' if \code{x} is tall (i.e., \code{nrow(x)}>>\code{ncol(x)} and
##' \code{ncol(x)}<30. If \code{ncol(x)}>30 then \code{matsort} simply calls
##' `\code{t(apply(x,1,sort))}'.  \code{matorder} returns a permutation which
##' rearranges its first argument into ascending order, breaking ties by
##' further arguments.
##' @author Raoul Grasman
##' @keywords array arith
##' @examples
##' 
##' example(Unique)
##' 
##' @export matmax matmin matsort matorder
"matmax" <-
function (...)
{
    x = cbind(...)
    if(!is.numeric(x))
        stop("Input should by numeric.")
    if (!is.matrix(drop(x)))
        x = t(x)
    x[1:nrow(x) + nrow(x) * (max.col(x) - 1)]
}
