dmvnorm <-
function (x, mu, Sigma, log = FALSE) {
    if (!is.matrix(x))
        x <- rbind(x)
    p <- length(mu)
    if (p == 1) {
        dnorm(x, mu, sqrt(Sigma), log = log)
    } else {
        t1 <- length(mu) == length(Sigma)
        t2 <- all(abs(Sigma[lower.tri(Sigma)]) < sqrt(.Machine$double.eps))
        if (t1 || t2) {
            if (!t1)
                Sigma <- diag(Sigma)
            nx <- nrow(x)
            ff <- rowSums(dnorm(x, rep(mu, each = nx), 
                sd = rep(sqrt(Sigma), each = nx), log = TRUE))
            if (log) ff else exp(ff)
        } else {
            ed <- eigen(Sigma, symmetric = TRUE)
            ev <- ed$values
            evec <- ed$vectors
            if (!all(ev >= -1e-06 * abs(ev[1]))) 
                stop("'Sigma' is not positive definite")
            ss <- x - rep(mu, each = nrow(x))
            inv.Sigma <- evec %*% (t(evec) / ev)
            quad <- 0.5 * rowSums((ss %*% inv.Sigma) * ss)
            fact <- - 0.5 * (p * log(2 * pi) + sum(log(ev)))
            if (log)
                as.vector(fact - quad)
            else
                as.vector(exp(fact - quad))
        }
    }
}
