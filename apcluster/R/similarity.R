negDistMat <- function(x, sel=NA, r=1, method="euclidean", p=2)
{
    if (r <= 0)
        stop("'r' must be positive")

    if (missing(x))
        return(function(x, sel=NA) negDistMat(x, sel, r=r, method=method, p=p))

    if (is.data.frame(x))
        dm <- as.matrix(simpleDist(x[, sapply(x, is.numeric)], sel,
                                   method=method, p=p))
    else
        dm <- as.matrix(simpleDist(x, sel, method=method, p=p))

    if (r != 1)
        -dm^r
    else
        -dm
}

expSimMat <- function(x, sel=NA, r=2, w=1, method="euclidean", p=2)
{
    if (r <= 0)
        stop("'r' must be positive")

    if (missing(x))
        return(function(x, sel=NA) expSimMat(x, sel, r=r, w=w,
                                             method=method, p=p))
    if (is.data.frame(x))
        dm <- as.matrix(simpleDist(x[,sapply(x, is.numeric)], sel,
                                   method=method, p=p))
    else
        dm <- as.matrix(simpleDist(x, sel, method=method, p=p))

    exp(-(dm / w)^r)
}

linSimMat <- function(x, sel=NA, w=1, method="euclidean", p=2)
{
    if (w <= 0)
        stop("'w' must be positive")

    if (missing(x))
        return(function(x, sel=NA) linSimMat(x, sel, w=w,
                                             method=method, p=p))
    if (is.data.frame(x))
        dm <- as.matrix(simpleDist(x[,sapply(x,is.numeric)], sel,
                                   method=method, p=p))
    else
        dm <- as.matrix(simpleDist(x, sel, method=method, p=p))

    pmax(1 - dm / w, 0)
}

corSimMat <- function(x, sel=NA, r=1, signed=TRUE, method="pearson")
{
    if (missing(x))
        return(function(x, sel=NA) corSimMat(x, sel, r=r, signed=signed,
                                             method=method))

    if (is.vector(x) || (is.list(x) && !is.data.frame(x)))
        stop("no correlation for vector or list")

    if (is.data.frame(x))
        x <- as.matrix(x[, sapply(x, is.numeric)])
    else
        x <- as.matrix(x)

    N  <- nrow(x)

    # if rownames available they are assigned by cor
    if (length(sel) == 1 && is.na(sel))
    {
        mat <- cor(x=t(x), method=method)

        if (length(rownames(x)) == 0)
            dimnames(mat) <- list(seq_len(N), seq_len(N))
    }
    else if (is.numeric(sel) && length(sel) > 0)
    {
        mat <- cor(x=t(x), y=t(x[sel, ]), method=method)

        if (length(rownames(x)) == 0)
            dimnames(mat) <- list(seq_len(N), sel)
    }
    else
        stop("invalid 'sel' argument")

    if (signed)
    {
        if (r != 1)
            mat <- sign(mat) * abs(mat)^r
    }
    else
    {
        if (r == 1)
            mat <- abs(mat)
        else
            mat <- abs(mat)^r
    }

    mat
}

linKernel <- function(x, sel=NA, normalize=FALSE)
{
    if (missing(x))
        return(function(x, sel=NA) linKernel(x, sel, normalize=normalize))

    if (is.data.frame(x))
        x <- as.matrix(x[, sapply(x, is.numeric)])
    else
        x <- as.matrix(x)

    N  <- nrow(x)

    if (!is.double(x)) storage.mode(x) <- "double"

    if (length(sel) == 1 && is.na(sel))
    {
        mat <- tcrossprod(x)

        if (normalize)
        {
            di <- 1 / sqrt(diag(mat))
            di[which(is.infinite(di))] <- 0

            mat <- mat * (di %o% di)
        }

        if (length(rownames(x)) > 0)
            dimnames(mat) <- list(rownames(x), rownames(x))
        else
            dimnames(mat) <- list(seq_len(N), seq_len(N))
    }
    else if (is.numeric(sel) && length(sel) > 0)
    {
        mat <- tcrossprod(x, x[sel, , drop=FALSE])

        if (normalize)
        {
            di <- 1 / sqrt(sapply(1:nrow(x), function(i) x[i,] %*% x[i,]))
            di[which(is.infinite(di))] <- 0

            mat <- mat * (di %o% di[sel])
        }

        if (length(rownames(x)) > 0)
            dimnames(mat) <- list(rownames(x), rownames(x)[sel])
        else
            dimnames(mat) <- list(seq_len(N), sel)
    }
    else
        stop("invalid 'sel' argument")

    mat
}
