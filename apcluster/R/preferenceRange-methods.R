preferenceRange.matrix <- function(s, exact=FALSE)
{
    if (ncol(s) != nrow(s))
        stop("'s' must be a square matrix")

    diag(s) <- 0

    if (length(which(s == -Inf)) > 0)
        warning("similarity matrix 's' contains -Inf similarities; ",
                "lower bound may not correspond to 1 or 2 clusters")

    .Call("preferenceRangeC", s, as.logical(exact)[1], PACKAGE="apcluster")
 }

setMethod("preferenceRange", signature(s="matrix"), preferenceRange.matrix)


preferenceRange.dgTMatrix <- function(s, exact=FALSE)
{
    if (ncol(s) != nrow(s))
        stop("'s' must be a square matrix")

    ## remove diagonal elements and -Inf from s
    remElem <- which(s@i == s@j | s@x == -Inf)

    if (length(remElem) > 0)
    {
        s@i <- s@i[-remElem]
        s@j <- s@j[-remElem]
        s@x <- s@x[-remElem]
    }

    .Call("preferenceRangeSparseC", s@i, s@j, s@x, nrow(s),
          as.logical(exact)[1], PACKAGE="apcluster")
 }

setMethod("preferenceRange", signature(s="dgTMatrix"),
          preferenceRange.dgTMatrix)


preferenceRange.otherSparse <- function(s, exact=FALSE)
{
    s <- try(as(as(s, "TsparseMatrix"), "dgTMatrix"))

    if (class(s) == "try-error")
        stop("cannot cast 's' (class '", class(s), "') to class 'dgTMatrix'")

    preferenceRange.dgTMatrix(s=s, exact=exact)
}

setMethod("preferenceRange", signature(s="sparseMatrix"),
          preferenceRange.otherSparse)


preferenceRange.otherDense <- function(s, exact=FALSE)
{
    s <- try(as(s, "matrix"))

    if (class(s) == "try-error")
        stop("cannot cast 's' (class '", class(s), "') to class 'matrix'")

    preferenceRange.matrix(s=s, exact=exact)
}

setMethod("preferenceRange", signature(s="Matrix"),
          preferenceRange.otherDense)
