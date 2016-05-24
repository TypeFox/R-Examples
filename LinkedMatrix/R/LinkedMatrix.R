#' @include ColumnLinkedMatrix.R RowLinkedMatrix.R
NULL


#' Initializes either a
#' \code{\link[=ColumnLinkedMatrix-class]{ColumnLinkedMatrix}} or
#' \code{\link[=RowLinkedMatrix-class]{RowLinkedMatrix}} instance of certain
#' dimensions with a configurable number and type of nodes.
#'
#' @param nrow The number of rows.
#' @param ncol The number of columns.
#' @param nNodes The number of nodes.
#' @param linkedBy Whether the matrix is linked by \code{rows} or
#'   \code{columns}.
#' @param nodeInitializer The name of a function or a function with four
#'   parameters \code{nodeIndex}, \code{ncol}, \code{nrow}, and \code{...} that
#'   initializes each node by returning a matrix-like object Pre-defined node
#'   initializers include \code{matrixNodeInitializer} to initialize matrices
#'   and \code{ffNodeInitializer} to initialize \code{ff} objects.
#' @param ... Additional arguments passed into \code{nodeInitializer}.
#' @export
LinkedMatrix <- function(nrow, ncol, nNodes, linkedBy, nodeInitializer, ...) {
    class <- ifelse(linkedBy == "columns", "ColumnLinkedMatrix", "RowLinkedMatrix")
    # Look for an internal function first
    ex <- try(nodeInitializer <- get(nodeInitializer), silent = TRUE)
    if (class(ex) == "try-error") {
        nodeInitializer <- match.fun(nodeInitializer)
    }
    linkedMatrix <- new(class)
    ranges <- chunkRanges(ifelse(class == "ColumnLinkedMatrix", ncol, nrow), nNodes)
    for (i in seq_len(nNodes)) {
        if (class == "RowLinkedMatrix") {
            n <- ranges[2, i] - ranges[1, i] + 1
            p <- ncol
        } else {
            n <- nrow
            p <- ranges[2, i] - ranges[1, i] + 1
        }
        linkedMatrix[[i]] <- nodeInitializer(nodeIndex = i, nrow = n, ncol = p, ...)
    }
    return(linkedMatrix)
}


matrixNodeInitializer <- function(nodeIndex, nrow, ncol, ...) {
    matrix(nrow = nrow, ncol = ncol, ...)
}


ffNodeInitializer <- function(nodeIndex, nrow, ncol, vmode, ...) {
    if (!requireNamespace("ff", quietly = TRUE)) {
        stop("The ff package is needed for this function to work. Please install it.", call. = FALSE)
    }
    ff::ff(dim = c(nrow, ncol), vmode = vmode, ...)
}


show <- function(object) {
    d <- dim(object)
    cat(d[1], "x", d[2], "linked matrix of class", class(object), "\n")
}


#' @export
length.LinkedMatrix <- function(x) {
    prod(dim(x))
}


#' Converts a \code{\link[=LinkedMatrix-class]{LinkedMatrix}} instance to a
#' \code{matrix} (if small enough).
#'
#' @param x Either a \code{\link[=ColumnLinkedMatrix-class]{ColumnLinkedMatrix}}
#'   or a \code{\link[=RowLinkedMatrix-class]{RowLinkedMatrix}} object.
#' @param ... Additional arguments (unused).
#' @export
as.matrix.LinkedMatrix <- function(x, ...) {
    x[, , drop = FALSE]
}


#' Returns the number of nodes.
#'
#' @param x Either a \code{\link[=ColumnLinkedMatrix-class]{ColumnLinkedMatrix}}
#'   or a \code{\link[=RowLinkedMatrix-class]{RowLinkedMatrix}} object.
#' @export
nNodes <- function(x) {
    length(slot(x, ".Data"))
}


#' Returns the column or row indexes at which each node starts and ends.
#'
#' @param x Either a \code{\link[=ColumnLinkedMatrix-class]{ColumnLinkedMatrix}}
#'   or a \code{\link[=RowLinkedMatrix-class]{RowLinkedMatrix}} object.
#' @return A matrix.
#' @export
nodes <- function(x) {
    UseMethod("nodes")
}


#' Maps each column or row index of a linked matrix to the column or row index
#' of its corresponding node.
#'
#' @param x Either a \code{\link[=ColumnLinkedMatrix-class]{ColumnLinkedMatrix}}
#'   or a \code{\link[=RowLinkedMatrix-class]{RowLinkedMatrix}} object.
#' @return A matrix.
#' @export
index <- function(x) {
    UseMethod("index")
}


#' An abstract S4 class union of
#' \code{\link[=ColumnLinkedMatrix-class]{ColumnLinkedMatrix}} and
#' \code{\link[=RowLinkedMatrix-class]{RowLinkedMatrix}}.
#'
#' This class is a class union and can therefore not be initialized. It can be
#' used to check whether an object is either of type
#' \code{\link[=ColumnLinkedMatrix-class]{ColumnLinkedMatrix}} or of type
#' \code{\link[=RowLinkedMatrix-class]{RowLinkedMatrix}} using \code{is(x,
#' "LinkedMatrix")}, and to assign methods for both
#' \code{\link[=ColumnLinkedMatrix-class]{ColumnLinkedMatrix}} and
#' \code{\link[=RowLinkedMatrix-class]{RowLinkedMatrix}} classes, e.g.
#' \code{\link[=show,LinkedMatrix-method]{show}}.
#'
#' @name LinkedMatrix-class
#' @docType class
#' @seealso \code{\link[=ColumnLinkedMatrix-class]{ColumnLinkedMatrix}} or
#'   \code{\link[=RowLinkedMatrix-class]{RowLinkedMatrix}} for implementations
#'   of column-linked matrices or row-linked matrices, respectively.
#' @exportClass LinkedMatrix
setClassUnion("LinkedMatrix", c("ColumnLinkedMatrix", "RowLinkedMatrix"))


#' Show a \code{\link[=LinkedMatrix-class]{LinkedMatrix}} object.
#'
#' This method is run when a \code{\link[=LinkedMatrix-class]{LinkedMatrix}}
#' object is printed.
#'
#' @param object Either a
#'   \code{\link[=ColumnLinkedMatrix-class]{ColumnLinkedMatrix}} or a
#'   \code{\link[=RowLinkedMatrix-class]{RowLinkedMatrix}} object.
#' @export
setMethod("show", signature(object = "LinkedMatrix"), show)
