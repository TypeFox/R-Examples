## matrices.R -- Part of the sparseHessianFD package
## Copyright (C) 2013-2016 Michael Braun


#' @name Matrix.to.Coord
#' @title Row and column indices from sparse matrix.
#' @description Utility function to extract row and column indices of
#' the non-zero
#' elements of a sparse matrix.
#' @param M A matrix that is a subclass of sparseMatrix, as defined in
#' the Matrix package.
#' @param index1 TRUE if the index of the first element should be 1,
#' and FALSE if 0.
#' @return A list with two named elements.
#' \describe{
#' \item{rows}{Integer vector containing row indices of non-zero elements}
#' \item{cols}{Integer vector containing column indices of non-zero elements}
#' }
#' @details A wrapper to \link{Matrix.to.Pointers} for \code{order='triplet'}
#' and \code{values=FALSE}, for extracting the row and column indices
#' of a sparsity pattern from a matrix that has that same pattern.
#' @examples
#' M1 <- as(kronecker(diag(3), matrix(TRUE,2,2)),"sparseMatrix")
#' C <- Matrix.to.Coord(M1)
#' M2 <- Matrix::sparseMatrix(i=C$rows, j=C$cols)
#' all.equal(M1,M2)
#' @export
Matrix.to.Coord <- function(M, index1=TRUE) {
    res <- Matrix.to.Pointers(M, values=FALSE, order="triplet", index1=index1)
    list(rows=res$rows, cols=res$cols)
}


#' @name Coord.to.Pointers
#' @title Convert a matrix defined by row and column indices to one
#' defined by a row- or column-oriented compression scheme.
#' @description Returns indices and pointers that define a sparse
#' Hessian in compressed format.  Inputs are the row and column indices.
#' @param rows,cols row and column indices of non-zero elements
#' @param dims 2-element vector for number of rows and columns.
#' @param order Determines the indexing/compression scheme for the
#' output matrix.  Use "triplet" to get row and column indices.
#' Defaults to the same class as M.
#' @param triangle Is input intended to be a triangular (TRUE) or full
#' (FALSE) matrix. See details for how this argument is interpreted
#' for different values of \code{order}.
#' @param lower If \code{triangular} is true, this argument identifies
#' the input matrix as lower- or upper-triangular.  This argument is
#' ignored if \code{triangle} is FALSE.
#' @param symmetric If TRUE, and matrix is triangular, then the matrix
#' will be treated as symmetric, with only the triangular elements
#' provided.  If matrix is neither triangular nor symmetric, then
#' symmetric=TRUE will probably trigger an error.
#' @param index1 TRUE if using 1-based indexing.  FALSE for 0-based indexing.
#' @details \code{triangle} and \code{order} have the following interpretation:
#' \describe{
#' \item{triangle=TRUE}{Input \code{rows} and {cols} represent lower
#' or upper triangle of a matrix. If \code{order="symmetric"}, then
#' the output list will be for a full, symmetric matrix. Otherwise,
#' the output list will be for only the lower or upper triangle.  Any
#' elements outside of the specified triangle will trigger an error.}
#' \item{triangle=FALSE}{Input \code{rows} and {cols} represent a full
#' matrix. If that matrix is not symmetric, then
#' \code{order=="symmetric"} will trigger an error.}
#' If \code{symmetric=FALSE} and \code{order='triplet'}, the output
#' list should contain the same row and column indices as the input
#' list.}
#' @return A list.  See Matrix.to.Pointers (no values are included in
#' return list).
#' @seealso Matrix.to.Pointers
#' @export
Coord.to.Pointers <- function(rows, cols, dims,
                              triangle=TRUE, lower=TRUE,
                              symmetric=FALSE,
                              order=c("column", "row", "triplet"),
                              index1=TRUE) {

    stopifnot(is.logical(triangle),
              is.logical(index1),
              is.logical(lower),
              length(rows) == length(cols)
              )

    if (triangle) {
        if (lower) {
            stopifnot(all(rows >= cols))
        } else {
            stopifnot(all(rows <= cols))
        }
    }

    R <- as(sparseMatrix(i=rows, j=cols, dims=dims,
                         index1=index1, giveCsparse=TRUE), "nMatrix")
    if (triangle & symmetric) {
        C <- as(sparseMatrix(i=cols, j=rows, dims=dims, index1=index1), "nMatrix")
        A <- as(R + C, "ngCMatrix") ## symmetric , but stored as general CSC sparse
    } else {
        A <- R
    }
    Matrix.to.Pointers(A, symmetric, values=FALSE, order=order, index1=index1)
}


#' @name Matrix.to.Pointers
#' @title Extract row and column indices, pointers and values from a
#' sparse matrix.
#' @description Returns a list of row indices, column indices,
#' pointers, and/or values of a sparse Hessian.
#' @param M A sparse Matrix, as defined in the Matrix package.
#' @param as.symmetric Defaults to isSymmetric(M).  If M is symmetric,
#' and as.symmetric is FALSE, then index/pointer elements in the
#' output list will be labeled according to order. If M is not
#' symmetric, and as.symmetric is TRUE, then an error will be triggered.
#' @param values  If TRUE, values are returned in list as 'x'.
#' Defaults to TRUE for numeric and logical matrices, and FALSE for
#' pattern matrices.  If M is a pattern matrix, values=TRUE will
#' trigger a warning.
#' @param order Determines the indexing/compression scheme for the
#' output matrix.  Use 'triplet" to get row and column indices.
#' Defaults to the same class as M.
#' @param index1 TRUE (default) if return indices and pointers should
#' use 1-based indexing.  FALSE for 0-based indexing.
#' @details This function is included primarily for debugging
#' purposes.  It is used internally, but would not ordinarily be
#' called by an end user.
#' @return A list with the following elements. If order=='row',
#' \describe{
#' \item{jCol}{ Integer vector containing column indices of non-zero elements}
#' \item{ipntr}{ Integer vector containing pointers to elements of
#' jCol at which the next row begins.}
#' }
#' If order=='column'
#'  \describe{
#' \item{iRow}{ Integer vector containing row indices of non-zero elements}
#' \item{jpntr}{ Integer vector containing pointers to elements of
#' iRow at which the next column begins.}
#' }
#' If order=='triplet'
#' \describe{
#' \item{rows}{Row indices of non-zero elements}
#' \item{cols}{Column indices of non-zero elements}
#' }
#' If as.symmetric is TRUE, then the row/column orientation does not matter.
#'  \describe{
#' \item{idx}{ Integer vector containing indices of non-zero elements}
#' \item{pntr}{ Integer vector containing pointers to elements of idx
#' at which the next row or column begins.}
#' }
#' If values=TRUE, the return list includes x, the values of the
#' non-zero elements.
#' The 'class' element is the name of the sparse matrix class to which
#' the output corresponds (identifies numeric type, pattern, and
#' indexing/compression scheme).
#' @seealso Matrix.to.Coord
#' @export
Matrix.to.Pointers <- function(M,
                               as.symmetric = Matrix::isSymmetric(M),
                               values=!is(M,"nMatrix"),
                               order=NULL,
                               index1=TRUE) {

    stopifnot(is.logical(as.symmetric),
              is.logical(values),
              is.logical(index1))
    stopifnot((is.null(order) | order %in% c("triplet","column","row")))


    Y <- as(as(M,"sparseMatrix"), "generalMatrix")
    if (!values) {
        Y <- as(Y, "nMatrix")
    } else {
        if (is(Y, "nMatrix")) {
            warning("M is a pattern matrix.  values=TRUE is ignored.")
        }
    }

    comp <- NULL
    if (is.null(order)) {
        if (is(Y, "TsparseMatrix")) comp <- "triplet"
        if (is(Y, "CsparseMatrix")) comp <- "column"
        if (is(Y, "RsparseMatrix")) comp <- "row"
    } else {
        comp <- order
        if (comp=="triplet") Y <- as(Y, "TsparseMatrix")
        if (comp=="column") Y <- as(Y, "CsparseMatrix")
        if (comp=="row") Y <- as(Y, "RsparseMatrix")
    }

    res <- vector("list", 3+values)
    if (comp=="triplet") {
        names(res)[1:2] <- c("rows","cols")
        res$rows <- Y@i + index1
        res$cols <- Y@j + index1
    } else {
        if (as.symmetric) {
            stopifnot(Matrix::isSymmetric(as(Y,"CsparseMatrix")))
            names(res)[1:2] <- c("idx","pntr")
            if (is(Y, "CsparseMatrix")) {
                res$idx <- Y@i + index1
            }
            if (is(Y, "RsparseMatrix")) {
                res$idx <- Y@j + index1
            }
            res$pntr <- Y@p + index1
        } else {
            if (comp=="column") {
                names(res)[1:2] <- c("iRow","jpntr")
                res$iRow <- Y@i + index1
                res$jpntr <- Y@p + index1
            } else {
                names(res)[1:2] <- c("jCol","ipntr")
                res$jCol <- Y@j + index1
                res$ipntr <- Y@p + index1
            }
        }
    }

    if (values) {
        names(res)[3] <- "x"
        res$x <- Y@x
    }
    names(res)[3+values] <- "class"
    res$class <- class(Y)

    return(res)
}

