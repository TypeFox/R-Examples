#### Cluster - Internal Utilities
#### ============================ (new by Martin Maechler)

## This was size(); seems slightly useful in general
sizeDiss <- function(d)
{
    ## find 'n' for d == dissimilarity-like(<n obs.>), i.e. length(d)= n(n-1)/2
    discr <- 1 + 8 * length(d)
    sqrtdiscr <- round(sqrt(discr))
    if(sqrtdiscr^2 == discr) (1 + sqrtdiscr)/2 else NA
}

##' Return indices to *permute*  "dissimilarity" / "dist" entries for C (ex-Fortran) code setup
##'
##' Currently always used as:
##'   n <- attr(x, "Size")
##'   dv <- x[lower.to.upper.tri.inds(n)]
##' -->> FIXME: eventually do the above directly in C
##' @param n "Size" = number of objects, underlying the dist/dissimilarity
##' used in  ./agnes.q, ./clara.q,  ./diana.q  und ./pam.q :
##' *somewhat* related to Matrix:::indTri()
lower.to.upper.tri.inds <- function(n)
{
    n1 <- as.integer(n - 1)
    if(n1 < 1) stop("'n' must be >= 2")
    else if(n1 == 1) 1L
    else rep(seq_len(n1), seq_len(n1)) +
        c(0L, unlist(lapply(2:n1, function(k) cumsum(c(0L, (n - 2L):(n - k))))))
}

upper.to.lower.tri.inds <- function(n)
{
    if((n2 <- as.integer(n - 2L)) < 0) stop("'n' must be >= 2")
    rep(1L + cumsum(0:n2), (n - 1):1) +
	unlist(lapply(0:n2, function(k) cumsum(k:n2)))
}


meanabsdev <- function(y) mean(abs(y - mean(y, na.rm = TRUE)), na.rm = TRUE)
