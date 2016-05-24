## deprecated.R -- Part of the sparseHessianFD package
## Copyright (C) 2013-2016 Michael Braun

## Place to hold functions that will not be maintained in the future


#' @name deprecated
#' @aliases Sym.CSC.to.Matrix Coord.to.Sym.Pattern.Matrix
#' @title Deprecated functions
#' @description These functions were in earlier versions, but will no
#' longer be maintained, and are not even guaranteed to work now.
NULL

## #' @title Sym.CSC.to.Matrix
#' @description Build sparse matrix from data in CSC (column
#' compressed) format.
#' @param H a list containing Hessian data.  See details.
#' @param nvars the number of rows (and columns) in the matrix.
#' @return An object of Matrix class.
#' @details Use Matrix::sparseMatrix instead of Sym.CSC.to.Matrix.
#' @rdname deprecated
#' @export
Sym.CSC.to.Matrix <- function(H,nvars) {

    .Deprecated("Matrix::spMatrix")

  M <- new("dsCMatrix", i = H$indrow, p = H$jpntr, x = H$vals, Dim=c(nvars, nvars),uplo="L")
  return(M)
}



#' @name Coord.to.Sym.Pattern.Matrix
#' @inheritParams Sym.CSC.to.Matrix
#' @details
#' Use Coord.to.Pattern.Matrix with symmetric=TRUE instead of Coord.to.Sym.Pattern.Matrix.
#' @rdname deprecated
#' @export
Coord.to.Sym.Pattern.Matrix <- function(H, nvars) {

    .Deprecated("Matrix::sparseMatrix")

## coords are for lower triangle, but coerces to symmetric pattern matrix
## H is a list with two integer vectors:  iRow and jCol



    res <- new("nsTMatrix",i=as.integer(H$iRow-1), j=as.integer(H$jCol-1),
               Dim=c(as.integer(nvars), as.integer(nvars)),uplo="L")

    return(res)

}

## #' @name Coord.to.Pattern.Matrix
## #' @aliases Coord.to.Pattern.Matrix
#' @title Pattern matrix from row and column indices.
#' @description Converts row and column indices to a pattern Matrix
#' object of Matrix class
#' @param rows,cols row and column indices of non-zero elements
#' @param dims 2-element vector for number of rows and columns in
#' matrix
#' @param compressed If TRUE, returns a matrix is compressed column (default=TRUE)
#' @param symmetric If TRUE, matrix will be symmetric, and only the
#' lower triangular elements need to be provided (default=FALSE)
#' @param index1 TRUE if input row and col use 1-based indexing, and FALSE for 0-based indexing.
#' @return A sparse pattern matrix
#' @details This function is useful to prototype a sparsity pattern.
#' No assumptions are made about symmetry.
#' @rdname deprecated
#' @export
Coord.to.Pattern.Matrix <- function(rows, cols, dims, compressed=TRUE,
                                    symmetric=FALSE, index1=TRUE) {

    res <- sparseMatrix(i=as.integer(rows),
                        j=as.integer(cols),
                        dims=dims,
                        giveCsparse=compressed,
                        symmetric=symmetric,
                        index1=index1)
    return(res)
}




## #' @name new.sparse.hessian.obj
#' @title Deprecated constructor
#' @param x variable vector for initialization
#' @param fn R function that returns function value
#' @param gr R function that returns the gradient of the function
#' @param hs list of two vectors:  row and column indices of non-zero
#' elements of lower triangle of Hessian.  See details.
#' @param fd.method If TRUE, use direct method for computatation.  Otherwise, use
#' indirect/substitution method.  See references.
#' @param eps The perturbation amount for finite differencing of the
#' gradient to compute the Hessian. Defaults to
#' sqrt(.Machine$double.eps).
#' @param ... Other parameters to be passed to fn and gr.
#' @details hs is a list of two elements:
#' \describe{
#' \item{iRow}{ Integer vector of row indices of non-zero elements in
#' lower triangle of Hessian.}
#' \item{jCol}{ Integer vector of column indices of non-zero elements in
#' lower triangle of Hessian.}
#' }
#' @rdname deprecated
#' @export
new.sparse.hessian.obj <- function(x, fn, gr, hs, fd.method=0L, eps=sqrt(.Machine$double.eps),...) {

    .Deprecated("sparseHessianFD")
    if (is.null(hs))
        stop("sparseHessianFD: you must supply structure of the Hessian.")
    if (!is.list(hs))
        stop ("sparseHessianFD:  hs must be a list")
    if (!all.equal(names(hs), c("rows","cols"))) {
        if (all.equal(names(hs), c("iRow","jCol"))) {
            names(hs) <- c("rows","cols")
        }
    }
    if (!all.equal(names(hs),c("rows","cols"))) {
        stop ("sparseHessianFD:  Names of hs must be either (\"iRow, jCol\") or (\"rows, cols\")")
    }
    direct <- as.logical(fd.method)
    sparseHessianFD(x, fn=fn, gr=gr, rows=hs$rows, cols=hs$cols,
                    eps, index1=TRUE, ...)
}


#' @name sparseHessianFD.new
#' @title Create and initialize a new sparseHessianFD object
#' @details This function is deprecated.  Use \code{sparseHessianFD} instead.
#' @description This function is deprecated.  Use \code{sparseHessianFD} instead.
#' @param direct If TRUE, use direct method for computatation.  Otherwise, use
#' indirect/substitution method.  See references.
#' @return An object of class sparseHessianFD
#' @rdname deprecated
#' @export
sparseHessianFD.new <- function(x, fn, gr, rows, cols, direct=NULL,
                            eps=sqrt(.Machine$double.eps), ...) {

    .Deprecated("sparseHessianFD")
    sparseHessianFD(x, fn=fn, gr=gr,
                           rows=rows, cols=cols, delta=eps,
                           index1=TRUE, ...)
}


