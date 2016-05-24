# This file is intented to be called by calex_1d.R.  It specifies
# E.theta() for the 1-d case.

E.theta.1d <- function(D2 = NULL, H1 = NULL, x1 = NULL, x2 = NULL, phi, give.mean = TRUE)
{
    if (give.mean) {
        m_theta <- phi$theta.apriori$mean
        return(H1(D1.fun(D2, t.vec = m_theta)))
    }
    else {
        out <- matrix(0, 3, 3)
        rownames(out) <- c("const","x","A")
        colnames(out) <- c("const","x","A")
        return(out)
    }
}
# Edash.theta.1d is just the same as the toy version (because we don't
# need the variance bits):

Edash.theta.1d <- Edash.theta.toy 

