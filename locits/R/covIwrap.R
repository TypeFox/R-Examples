covIwrap <-
function (S, m, n, ll, storewrap, P) 
{
    if (missing(storewrap)) 
        stop("covI wrap has to have storewrap, which stores intermediate results")
    if (!is.null(storewrap)) {
        if (any(S != storewrap$S)) 
            stop("New S supplied is not same as stored")
        the.store <- storewrap$the.store
    }
    else the.store <- NULL
    d <- abs(m - n)
    if (!is.null(the.store)) {
        ix1 <- the.store[, 1] == ll
        ix2 <- the.store[, 2] == d
        ix <- ix1 & ix2
        nfound <- sum(ix)
        if (nfound > 1) 
            stop("Confused: why are there more than one found?")
        if (nfound == 1) {
            l <- list(ans = the.store[ix, 3], storewrap = storewrap)
            return(l)
        }
    }
    ans <- covI(S, m, n, ll, ThePsiJ = P)
    if (is.null(the.store)) {
        the.store <- matrix(c(ll, d, ans), nrow = 1, ncol = 3)
        storewrap <- list(S = S)
    }
    else {
        the.store <- rbind(the.store, c(ll, d, ans))
    }
    storewrap$the.store <- the.store
    l <- list(ans = ans, storewrap = storewrap)
    l
}
