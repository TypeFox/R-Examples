##' Apply a Function to a Data Frame Split by Factors
##'
##' Simple wrapper of the 'by' function
##' @title Apply a Function to a Data Frame Split by Factors
##' @param x Data frame
##' @param INDICES Indices (vector or list of indices, vector of column names, or formula of column names)
##' @param FUN A function to be applied to data frame subsets of 'data'.
##' @param COLUMNS (Optional) subset of columns of x to work on
##' @param array if TRUE an array/matrix is always returned
##' @param ... Additional arguments to lower-level functions
##' @author Klaus K. Holst
##' @export
##' @examples
##' By(datasets::CO2,~Treatment+Type,colMeans,~conc)
##' By(datasets::CO2,~Treatment+Type,colMeans,~conc+uptake)
By <- function(x,INDICES,FUN,COLUMNS,array=FALSE,...) {
    if (inherits(INDICES,"formula")) {
        INDICES <- as.list(model.frame(INDICES,x))
    } else {
        if (is.character(INDICES) && length(INDICES)!=nrow(x)) {
            INDICES <- as.list(x[,INDICES,drop=FALSE])
        }
    }
    if (!missing(COLUMNS)) {
        if (inherits(COLUMNS,"formula")) {
            x <- model.frame(COLUMNS,x)
        } else {
            x <- x[,COLUMNS,drop=FALSE]
        }
    }
    a <- by(x, INDICES, FUN=FUN, ...)
    if (NCOL(x)==1 && !array) {
        ##DimElem <- length(a[rep(1,length(dim(a)))][[1]])
        a <- a[]
        attr(a,"call") <- NULL
        ##        a <- array(a,)
    }
    return(a)
}
