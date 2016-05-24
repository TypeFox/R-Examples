##
##  i n v l a p . R  Inverse Laplacian
##


invlap <- function(Fs, t1, t2, nnt, a = 6, ns = 20, nd = 19) {
    stopifnot(is.numeric(t1), length(t1) == 1, is.numeric(t2), length(t2) == 1,
              is.numeric(nnt), length(nnt) == 1)
    Fs <- match.fun(Fs)

    radt <- linspace(t1, t2, nnt)
    if (t1 == 0)  {
        radt <- radt[2:nnt]
        nnt <- nnt - 1
    }

    alfa <- beta <- numeric(ns+1+nd)
    for (n in 1:(ns+1+nd)) {
       alfa[n] <- a + (n-1) * pi * 1i
       beta[n] <- -exp(a) * (-1)^n
    }
    n <- 1:nd
    bdif <- rev(cumsum(gamma(nd+1) / gamma(nd+2-n) / gamma(n))) / 2^nd
    beta[(ns+2):(ns+1+nd)] <- beta[(ns+2):(ns+1+nd)] * bdif
    beta[1] = beta[1]/2

    ft <- numeric(nnt)
    for (kt in 1:nnt) {
       tt <- radt[kt]
       s <- alfa/tt
       bt <- beta/tt
       btF <- bt * Fs(s)
       ft[kt] <- sum(Re(btF))
    }

    return(list(x = radt, y = ft))
}
