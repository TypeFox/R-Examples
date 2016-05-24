##' Retrieve or set a list of array element values
##' 
##' \code{entry.value} retrieves or sets the values in an array \code{a} at the
##' positions indicated by the rows of a matrix \code{idx}.
##' 
##' 
##' @aliases entry.value entry.value<-
##' @param a An array.
##' @param idx Numerical matrix with the same number of columns as the number
##' of dimensions of \code{a}.  Each row indices a cell in \code{a} of which
##' the value is to be retrieved or set.
##' @param value An array of length \code{nrow(idx)}.
##' @return \code{entry.value(a,idx)} returns a vector of values at the
##' indicated cells.  \code{entry.value(a,idx) <- val} changes the indicated
##' cells of \code{a} to \code{val}.
##' @author Raoul Grasman
##' @keywords arith math array
##' @examples
##' 
##' a = array(1:(4^4),c(4,4,4,4))
##' entry.value(a,cbind(1:4,1:4,1:4,1:4))
##' entry.value(a,cbind(1:4,1:4,1:4,1:4)) <- 0
##' 
##' entry.value(a, as.matrix(expand.grid(1:4,1:4,1:4,1:4)))
##'      # same as `c(a[1:4,1:4,1:4,1:4])' which is same as `c(a)'
##'
##' @export entry.value entry.value<-
"entry.value" <-
function (a, idx) 
{
    if (!is.array(a)) 
        stop(paste("First argument `", deparse(substitute(a)), 
            "' should be an array.", sep = ""))
    if (!is.matrix(idx)) 
        stop(paste("Second argument `", substitute(idx), "' should be a matrix.", 
            sep = ""))
    n <- length(dim(a))
    if (n != ncol(idx)) 
        stop(paste("Number of columns in", deparse(substitute(idx)), 
            "is incompatible is dimension of", deparse(substitute(a))))
    a[(idx - 1) %*% c(1, cumprod(dim(a))[-n]) + 1]
}
"entry.value<-" <-
function (a, idx, value) 
{
    if (!is.array(a)) 
        stop(paste("First argument `", deparse(substitute(a)), 
            "' should be an array.", sep = ""))
    if (!is.matrix(idx)) 
        stop(paste("Second argument `", substitute(idx), "' should be a matrix.", 
            sep = ""))
    n <- length(dim(a))
    if (n != ncol(idx)) 
        stop(paste("Number of columns in", deparse(substitute(idx)), 
            "is incompatible is dimension of", deparse(substitute(a))))
    a[(idx - 1) %*% c(1, cumprod(dim(a))[-n]) + 1] <- value
    return(a)
}

