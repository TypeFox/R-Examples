dist1 <- 
function (x, X, p = 2, inf = FALSE)
{
    X <- if (is.null(dim(X)))
        matrix(X, ncol = length(x))
    else as.matrix(X)
    dif <- abs(X - matrix(as.numeric(x), nrow = nrow(X), ncol = ncol(X), byrow = 
        TRUE))
    if (inf)
        return(apply(dif, 1, max))
    tmp <- dif^p
    rowSums(tmp)
}
