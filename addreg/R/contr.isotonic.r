contr.isotonic <- function (n, perm, contrasts = TRUE, sparse = FALSE)
{
    if(is.numeric(n) && length(n) == 1L) {
        if (n > 1L) levels <- as.character(seq_len(n))
        else stop("not enough degrees of freedom to define contrasts")
    }
    else {
        levels <- as.character(n)
        n <- length(n)
    }
    if(!missing(perm)) {
        if(length(perm) != length(levels))
            stop("'perm' must have the same number of levels as 'n'")
        if(is.character(perm)) {
            if(any(sort(perm) != sort(levels)))
                stop("'perm' must be a permutation of the levels in 'n'")
            permn <- match(levels,perm)
        }
        else if(is.numeric(perm)) {
            if(any(sort(perm) != seq_len(n)))
                stop("'perm' must be a permutation of 'n'")
            permn <- perm
            perm <- levels[perm]
        }
    }
    else {
        permn <- seq_len(n)
        perm <- as.character(permn)
    }
    contr <- .Triang(perm)[permn,]
    if(contrasts) contr[,-1, drop = FALSE]
    else contr
}

.Triang <- function(nms) {
    n <- as.integer(length(nms))
    d <- c(n,n)
    dn <- list(nms, nms)
    z <- array(0,d,dn)
    z[row(z)>=col(z)] <- 1L
    z
}