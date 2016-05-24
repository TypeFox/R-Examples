# This is a modification of
#     File src/library/stats/R/dist.R
#     Part of the R package, http://www.R-project.org
# Changes:
#     added possibility to make sub-selection of columns (argument sel)
#     removed all functionality not necessary for package apcluster
#  This program is free software; you can redistribute it and/or modify

simpleDist <- function(x, sel, method="euclidean", p=2)
{
    ## account for possible spellings of euclid?an
    if(!is.na(pmatch(method, "euclidian")))
	method <- "euclidean"

    METHODS <- c("euclidean", "maximum",
		 "manhattan", "canberra", "binary", "minkowski")
    method <- pmatch(method, METHODS)
    if(is.na(method))
	stop("invalid distance method")
    if(method == -1)
	stop("ambiguous distance method")

    x <- as.matrix(x)
    N  <- nrow(x)

    if (!is.double(x)) storage.mode(x) <- "double"

    if (length(sel) == 1 && is.na(sel))
    {
        d <- .Call("CdistR", x, as.integer(NA), method, p)

        dm <- matrix(0, N, N)
        dm[row(dm) > col(dm)] <- d

        dm <- dm + t(dm)

        if (length(rownames(x)) > 0)
            dimnames(dm) <- list(rownames(x), rownames(x))
        else
            dimnames(dm) <- list(seq_len(N), seq_len(N))
    }
    else if (is.numeric(sel) && length(sel) > 0)
    {
        if (max(sel) > N || min(sel) < 1)
            stop("'sel' is no subset of '1:nrow(x)'")

        d <- .Call("CdistR", x, as.integer(sel - 1), method, p)

        dm <- matrix(d, N, length(sel))

        if (length(rownames(x)) > 0)
            dimnames(dm) <- list(rownames(x), rownames(x)[sel])
        else
            dimnames(dm) <- list(seq_len(N), sel)
    }
    else
        stop("invalid 'sel' argument")

    dm
}
