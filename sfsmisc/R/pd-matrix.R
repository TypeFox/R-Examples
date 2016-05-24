
## testing code is currently in
## /u/maechler/R/MM/MISC/posdefify.R

## TODO: probaby add 'rescale.kind = c("diag", "trace")'
## ----  they would differ only when some EV's were negative

### Higham's code by Ravi

posdefify <- function(m, method = c("someEVadd", "allEVadd"),
                      symmetric = TRUE,
                      eigen.m = eigen(m, symmetric= symmetric),
                      eps.ev = 1e-7)
{
    ## Purpose: From a matrix m, make a "close" positive definite one
    ## -------------------------------------------------------------------------
    ## Arguments: m: numeric matrix (n x n), usually symmetric
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 19 Dec 1997; 7 Jul 2004
    stopifnot(is.numeric(m) && is.matrix(m))
    method <- match.arg(method)
    n <- length(lam <- eigen.m $values)
    Eps <- eps.ev * abs(lam[1])# lam[1] is largest EV; "small" is *relative*
    ## lam[n] is the SMALLEST eigenvalue
    if(lam[n] < Eps) { # fix up small or negative values
        switch(method,
               "someEVadd" = lam[lam < Eps] <- Eps,
               "allEVadd" =  lam <- lam + Eps-lam[n]
               )
        Q <- eigen.m $vectors
        o.diag <- diag(m)# original one - for rescaling
        m <- Q %*% (lam * t(Q)) ## == Q %*% diag(lam) %*% t(Q)
        ## rescale to the original diagonal values
        ## D <- sqrt(o.diag/diag(m))
        ## where they are >= Eps :
        D <- sqrt(pmax(Eps,o.diag)/diag(m))
        m[] <- D * m * rep(D, each = n) ## == diag(D) %*% m %*% diag(D)
    }
    m
}
