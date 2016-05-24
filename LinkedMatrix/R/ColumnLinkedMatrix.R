subset.ColumnLinkedMatrix <- function(x, i, j, ..., drop) {
    nX <- nrow(x)
    pX <- ncol(x)
    if (missing(i)) {
        i <- 1:nX
    }
    if (missing(j)) {
        j <- 1:pX
    }
    if (class(i) == "logical") {
        i <- rep_len(i, nX)
        i <- which(i)
    } else if (class(i) == "character") {
        i <- sapply(i, function(name) {
            which(rownames(x) == name)
        }, USE.NAMES = FALSE)
    }
    if (class(j) == "logical") {
        j <- rep_len(j, pX)
        j <- which(j)
    } else if (class(j) == "character") {
        j <- sapply(j, function(name) {
            which(colnames(x) == name)
        }, USE.NAMES = FALSE)
    }
    n <- length(i)
    p <- length(j)
    originalOrder <- (1:p)[order(j)]
    sortedColumns <- sort(j)
    if (p > pX | n > nX) {
        stop("Either the number of columns or number of rows requested exceed the number of rows or columns in x, try dim(x)...")
    }
    Z <- matrix(nrow = n, ncol = p, NA)
    colnames(Z) <- colnames(x)[j]
    rownames(Z) <- rownames(x)[i]
    INDEX <- index(x)[sortedColumns, , drop = FALSE]
    whatChunks <- unique(INDEX[, 1])
    end <- 0
    for (k in whatChunks) {
        TMP <- matrix(data = INDEX[INDEX[, 1] == k, ], ncol = 3)
        ini <- end + 1
        end <- ini + nrow(TMP) - 1
        # Convert to matrix to support data frames
        Z[, ini:end] <- as.matrix(x[[k]][i, TMP[, 3], drop = FALSE])
    }
    if (length(originalOrder) > 1) {
        Z[] <- Z[, originalOrder]
    }
    if (drop == TRUE && (n == 1 || p == 1)) {
        # Revert drop.
        return(Z[, ])
    } else {
        return(Z)
    }
}


replace.ColumnLinkedMatrix <- function(x, i, j, ..., value) {
    if (missing(i)) {
        i <- 1:nrow(x)
    }
    if (missing(j)) {
        j <- 1:ncol(x)
    }
    Z <- matrix(nrow = length(i), ncol = length(j), data = value)
    # Retrieve nodes and index from ... to speed up sequential writes
    ellipsis <- list(...)
    if (is.null(ellipsis$nodes)) {
        nodes <- nodes(x)
    } else {
        nodes <- ellipsis$nodes
    }
    if (is.null(ellipsis$index)) {
        index <- index(x)
    } else {
        index <- ellipsis$index
    }
    for (k in 1:nrow(nodes)) {
        col_z <- (j >= nodes[k, 2]) & (j <= nodes[k, 3])
        colLocal <- index[j[col_z], 3]
        x[[k]][i, colLocal] <- Z[, col_z]
    }
    return(x)
}


#' @export
dim.ColumnLinkedMatrix <- function(x) {
    n <- nrow(x[[1]])
    p <- 0
    for (i in 1:nNodes(x)) {
        p <- p + ncol(x[[i]])
    }
    return(c(n, p))
}


# This function looks like an S3 method, but isn't one.
rownames.ColumnLinkedMatrix <- function(x) {
    out <- rownames(x[[1]])
    return(out)
}


# This function looks like an S3 method, but isn't one.
colnames.ColumnLinkedMatrix <- function(x) {
    out <- NULL
    if (!is.null(colnames(x[[1]]))) {
        p <- dim(x)[2]
        out <- rep("", p)
        nodes <- nodes(x)
        for (i in 1:nrow(nodes)) {
            out[(nodes[i, 2]:nodes[i, 3])] <- colnames(x[[i]])
        }
    }
    return(out)
}


#' @export
dimnames.ColumnLinkedMatrix <- function(x) {
    list(rownames.ColumnLinkedMatrix(x), colnames.ColumnLinkedMatrix(x))
}


# This function looks like an S3 method, but isn't one.
`rownames<-.ColumnLinkedMatrix` <- function(x, value) {
    for (i in 1:nNodes(x)) {
        rownames(x[[i]]) <- value
    }
    return(x)
}


# This function looks like an S3 method, but isn't one.
`colnames<-.ColumnLinkedMatrix` <- function(x, value) {
    nodes <- nodes(x)
    for (i in 1:nrow(nodes)) {
        colnames(x[[i]]) <- value[(nodes[i, 2]:nodes[i, 3])]
    }
    return(x)
}


#' @export
`dimnames<-.ColumnLinkedMatrix` <- function(x, value) {
    d <- dim(x)
    rownames <- value[[1]]
    colnames <- value[[2]]
    if (!is.list(value) || length(value) != 2 || !(is.null(rownames) || length(rownames) == d[1]) || !(is.null(colnames) ||
        length(colnames) == d[2])) {
        stop("invalid dimnames")
    }
    x <- `rownames<-.ColumnLinkedMatrix`(x, rownames)
    x <- `colnames<-.ColumnLinkedMatrix`(x, colnames)
    return(x)
}


#' Combine matrix-like objects by columns.
#'
#' Compared to the
#' \code{\link[=initialize,ColumnLinkedMatrix-method]{ColumnLinkedMatrix}}
#' constructor, nested \code{\link[=LinkedMatrix-class]{LinkedMatrix}} objects
#' that are passed via \code{...} will not be treated as matrix-like objects,
#' but their nodes will be extracted and merged with the new
#' \code{\link[=ColumnLinkedMatrix-class]{ColumnLinkedMatrix}} object for a
#' more compact representation.
#'
#' @param ... Matrix-like objects to be combined by columns.
#' @param deparse.level Currently unused, defaults to 0.
#' @export
cbind.ColumnLinkedMatrix <- function(..., deparse.level = 0) {
    dotdotdot <- list(...)
    nodes <- list()
    for (i in seq_len(length(dotdotdot))) {
        node <- dotdotdot[[i]]
        if (is(node, "LinkedMatrix")) {
            # Extract nodes from LinkedMatrix object
            nodes <- append(nodes, slot(node, ".Data"))
        } else {
            nodes <- append(nodes, node)
        }
    }
    do.call(ColumnLinkedMatrix, nodes)
}


#' Combine matrix-like objects by rows.
#'
#' @param ... Matrix-like objects to be combined by rows.
#' @param deparse.level Currently unused, defaults to 0.
#' @export
rbind.ColumnLinkedMatrix <- function(..., deparse.level = 1) {
    stop("rbind is currently undefined for ColumnLinkedMatrix")
}


#' @export
nodes.ColumnLinkedMatrix <- function(x) {
    n <- nNodes(x)
    OUT <- matrix(nrow = n, ncol = 3, NA)
    colnames(OUT) <- c("node", "col.ini", "col.end")
    end <- 0
    for (i in 1:n) {
        ini <- end + 1
        end <- ini + ncol(x[[i]]) - 1
        OUT[i, ] <- c(i, ini, end)
    }
    return(OUT)
}


#' @export
index.ColumnLinkedMatrix <- function(x) {
    nodes <- nodes(x)
    nColIndex <- nodes[nrow(nodes), 3]
    INDEX <- matrix(nrow = nColIndex, ncol = 3)
    colnames(INDEX) <- c("node", "col.global", "col.local")
    INDEX[, 2] <- 1:nColIndex
    end <- 0
    for (i in 1:nrow(nodes)) {
        nColChunk <- nodes[i, 3] - nodes[i, 2] + 1
        ini <- end + 1
        end <- ini + nColChunk - 1
        INDEX[ini:end, 1] <- i
        INDEX[ini:end, 3] <- 1:nColChunk
    }
    return(INDEX)
}


#' An S4 class to represent a column-linked
#' \code{\link[=LinkedMatrix-class]{LinkedMatrix}}.
#'
#' This class treats a list of matrix-like objects that are linked together by
#' columns and have the same number of rows similarly to a regular \code{matrix}
#' by implementing key methods such as \code{[} and \code{[<-} for extracting
#' and replacing matrix elements, \code{dim} to retrieve dimensions, and
#' \code{dimnames} and \code{dimnames<-} to retrieve and set dimnames. Each list
#' element is called a node and can be extracted or replaced using \code{[[} and
#' \code{[[<-}. A matrix-like object is one that has two dimensions and
#' implements at least \code{dim} and \code{[}.
#'
#' There are several ways to create an instance of this class: either by using
#' one of the constructors
#' \code{\link[=initialize,ColumnLinkedMatrix-method]{ColumnLinkedMatrix(...)}}
#' or
#' \code{\link[=initialize,ColumnLinkedMatrix-method]{new("ColumnLinkedMatrix",...)}},
#' or by using the more general \code{\link{LinkedMatrix}} function that
#' constructs objects of certain dimensions with a configurable number and type
#' of nodes.
#'
#' @export ColumnLinkedMatrix
#' @exportClass ColumnLinkedMatrix
ColumnLinkedMatrix <- setClass("ColumnLinkedMatrix", contains = "list")


#' Creates a new \code{\link[=ColumnLinkedMatrix-class]{ColumnLinkedMatrix}}
#' instance.
#'
#' This method is run when a
#' \code{\link[=ColumnLinkedMatrix-class]{ColumnLinkedMatrix}} object is created
#' using
#' \code{\link[=initialize,ColumnLinkedMatrix-method]{ColumnLinkedMatrix(...)}}
#' or
#' \code{\link[=initialize,ColumnLinkedMatrix-method]{new("ColumnLinkedMatrix",...)}}
#' and accepts a list of matrix-like objects as \code{...}. A matrix-like object
#' is one that has two dimensions and implements at least \code{dim} and
#' \code{[}. \code{\link[=LinkedMatrix-class]{LinkedMatrix}} objects can be
#' nested as long as they are conformable. Each object needs to have the same
#' number of rows to be linked together. If no matrix-like objects are given, a
#' single 1x1 node of type \code{matrix} filled with \code{NA} is returned.
#'
#' @inheritParams base::list
#' @param .Object The
#'   \code{\link[=ColumnLinkedMatrix-class]{ColumnLinkedMatrix}} instance to be
#'   initialized. This argument is passed in by R and can be ignored, but still
#'   needs to be documented.
#' @param ... A sequence of matrix-like objects of the same row-dimension.
#' @export
setMethod("initialize", signature(.Object = "ColumnLinkedMatrix"), function(.Object, ...) {
    nodes <- list(...)
    # Append at least one matrix
    if (length(nodes) == 0) {
        nodes[[1]] <- matrix()
    } else {
        # Detect non-matrix objects by checking dimensions
        if (any(sapply(nodes, function(x) length(dim(x)) != 2))) {
            stop("arguments need to be matrix-like")
        }
        # Detect matrices that do not match in dimensions
        if (length(unique(sapply(nodes, nrow))) != 1) {
            stop("arguments need the same number of rows")
        }
    }
    .Object <- callNextMethod(.Object, nodes)
    return(.Object)
})


#' Extract parts of a
#' \code{\link[=ColumnLinkedMatrix-class]{ColumnLinkedMatrix}}.
#'
#' This method is run when the \code{[]} operator is used on a
#' \code{\link[=ColumnLinkedMatrix-class]{ColumnLinkedMatrix}} object.
#'
#' @inheritParams base::`[`
#' @param j Column indices.
#' @param ... Additional arguments
#' @export
setMethod("[", signature(x = "ColumnLinkedMatrix"), subset.ColumnLinkedMatrix)


#' Replace parts of a
#' \code{\link[=ColumnLinkedMatrix-class]{ColumnLinkedMatrix}}.
#'
#' This method is run when the \code{[]} operator is used in an assignment on a
#' \code{\link[=ColumnLinkedMatrix-class]{ColumnLinkedMatrix}} object.
#'
#' @inheritParams base::`[<-`
#' @param j Column indices.
#' @param ... Additional arguments
#' @export
setReplaceMethod("[", signature(x = "ColumnLinkedMatrix"), replace.ColumnLinkedMatrix)
