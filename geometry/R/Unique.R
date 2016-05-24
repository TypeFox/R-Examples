##' Extract Unique Rows
##' 
##' \sQuote{Unique} returns a vector, data frame or array like 'x' but with
##' duplicate elements removed.
##' 
##' 
##' @param X Numerical matrix.
##' @param rows.are.sets If \sQuote{\code{TRUE}}, rows are treated as sets -
##' i.e., to define uniqueness, the order of the rows does not matter.
##' @return Matrix of the same number of columns as \code{x}, with the unique
##' rows in \code{x} sorted according to the columns of \code{x}. If
##' \code{rows.are.sets = TRUE} the rows are also sorted.
##' @note \sQuote{\code{Unique}} is (under circumstances) much quicker than the
##' more generic base function \sQuote{\code{unique}}.
##' @author Raoul Grasman
##' @keywords arith math array
##' @examples
##' 
##' # `Unique' is faster than `unique'
##' x = matrix(sample(1:(4*8),4*8),ncol=4)
##' y = x[sample(1:nrow(x),3000,TRUE), ]
##' gc(); system.time(unique(y))
##' gc(); system.time(Unique(y))
##' 
##' #
##' z = Unique(y)
##' x[matorder(x),]
##' z[matorder(z),]
##' 
##' @export
"Unique" <-
function (X, rows.are.sets = FALSE) 
{
    if (rows.are.sets) 
        X = matsort(X)
    X = X[matorder(X), ]
    dX = apply(X, 2, diff)
    uniq = c(TRUE, ((dX^2) %*% rep(1, ncol(dX))) > 0)
    X = X[uniq, ]
    return(X)
}
