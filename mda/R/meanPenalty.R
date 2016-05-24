meanPenalty <-
function (prior, cl) 
{
    n <- length(prior)
    Q <- matrix(0, n, n)
    ll <- unique(cl)
    for (lll in ll) {
        which <- cl == lll
        pp <- prior[which]
        npp <- length(pp)
        pp <- diag(npp) - outer(rep(1, npp), pp/sum(pp))
        Q[which, which] <- t(pp) %*% pp
    }
    attr(Q, "cl") <- cl
    attr(Q, "prior") <- prior
    Q
}

