sparseToFull <- function(s, fill=-Inf)
{
    .Deprecated("as.DenseSimilarityMatrix")

    as.DenseSimilarityMatrix(s=s, fill=fill)
}


as.SparseSimilarityMatrix.matrix <- function(s, lower=-Inf)
{
    if (!is(s, "matrix"))
    {
        s <- try(as(s, "matrix"))

        if (class(s) == "try-error")
            stop("cannot cast 's' (class '", class(s), "') to class 'matrix'")
    }

    if (nrow(s) != ncol(s) && ncol(s) != 3)
        stop("matrix 's' must be quadratic or have 3 columns")

    if (nrow(s) == ncol(s))
    {
        sel <- which(s > lower, arr.ind=TRUE)

        remElem <- which(sel[, 1] == sel[, 2])

        if (length(remElem) > 0)
            sel <- sel[-remElem, , drop=FALSE]

        if (nrow(sel) == 0)
            S <- new("dgTMatrix", Dim=dim(s))
        else
            S <- new("dgTMatrix", Dim=dim(s),
                     i=as.integer(sel[, 1] - 1),
                     j=as.integer(sel[, 2] - 1),
                     x=s[sel])
    }
    else
    {
        if (min(s[, 1:2]) <= 0)
            stop("indices in 's' must be >= 1")

        if (any(s[, 1:2] != floor(s[, 1:2])))
            stop("indices in 's' must be natural numbers")

        remElem <- which(s[, 1] == s[, 2] | s[, 3] <= lower)

        if (length(remElem) > 0)
            s <- s[-remElem, , drop=FALSE]

        if (nrow(s) == 0)
            S <- new("dgTMatrix", Dim=as.integer(c(0, 0)))
        else
        {
            N <- max(s[, 1:2])

            S <- new("dgTMatrix", Dim=as.integer(c(N, N)),
                     i=as.integer(s[, 1] - 1),
                     j=as.integer(s[, 2] - 1),
                     x=s[, 3])
        }
    }

    S
}

setMethod("as.SparseSimilarityMatrix", signature(s="matrix"),
          as.SparseSimilarityMatrix.matrix)

setMethod("as.SparseSimilarityMatrix", signature(s="Matrix"),
          as.SparseSimilarityMatrix.matrix)


as.SparseSimilarityMatrix.sparseMatrix <- function(s, lower=-Inf)
{
    if (nrow(s) != ncol(s))
        stop("argument 's' must be quadratic similarity matrix")

    if (!is(s, "dgTMatrix"))
    {
        s <- try(as(as(s, "TsparseMatrix"), "dgTMatrix"))

        if (class(s) == "try-error")
            stop("cannot cast 's' (class '", class(s),
                 "') to class 'dgTMatrix'")
    }

    remElem <- which(s@i == s@j | s@x <= lower)

    if (length(remElem) > 0)
    {
        s@i <- s@i[-remElem]
        s@j <- s@j[-remElem]
        s@x <- s@x[-remElem]
    }

    s
}

setMethod("as.SparseSimilarityMatrix", signature(s="sparseMatrix"),
          as.SparseSimilarityMatrix.sparseMatrix)



as.DenseSimilarityMatrix.matrix <- function(s, fill=-Inf)
{
    if (!is(s, "matrix"))
    {
        s <- try(as(s, "matrix"))

        if (class(s) == "try-error")
            stop("cannot cast 's' (class '", class(s), "') to class 'matrix'")
    }

    if (ncol(s) != 3)
        stop("'s' must be a matrix with 3 columns")

    if (min(s[, 1:2]) <= 0)
        stop("indices in 's' must be >= 1")

    if (any(s[, 1:2] != floor(s[, 1:2])))
        stop("indices in 's' must be natural numbers")

    N <- max(s[, 1:2])

    S <- matrix(fill, N, N)

    S[s[, 1] + N * (s[, 2] - 1)] <- s[, 3]

    S
}

setMethod("as.DenseSimilarityMatrix", signature(s="matrix"),
          as.DenseSimilarityMatrix.matrix)

setMethod("as.DenseSimilarityMatrix", signature(s="Matrix"),
          as.DenseSimilarityMatrix.matrix)


as.DenseSimilarityMatrix.sparseMatrix <- function(s, fill=-Inf)
{
    if (nrow(s) != ncol(s))
        stop("argument 's' must be quadratic similarity matrix")

    if (!is(s, "dgTMatrix"))
    {
        s <- try(as(as(s, "TsparseMatrix"), "dgTMatrix"))

        if (class(s) == "try-error")
            stop("cannot cast 's' (class '", class(s),
                 "') to class 'dgTMatrix'")
    }

    N <- nrow(s)

    S <- matrix(fill, N, N)

    S[(s@i + 1) + N * s@j] <- s@x

    S
}

setMethod("as.DenseSimilarityMatrix", signature(s="sparseMatrix"),
          as.DenseSimilarityMatrix.sparseMatrix)
