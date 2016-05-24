dmvt <-
function (x, mu, Sigma, df, log = FALSE) {
    if (!is.numeric(x))
        stop("'x' must be a numeric matrix or vector")
    if (!is.matrix(x))
        x <- rbind(x)
    p <- length(mu)
    if (!all(dim(Sigma) == c(p, p)) || ncol(x) != p) 
        stop("incompatible arguments")
    ed <- eigen(Sigma, symmetric = TRUE)
    ev <- ed$values
    if (!all(ev >= -1e-06 * abs(ev[1]))) 
        stop("'Sigma' is not positive definite")
    ss <- x - rep(mu, each = nrow(x))
    inv.Sigma <- ed$vectors %*% (t(ed$vectors) / ev)
    quad <- rowSums((ss %*% inv.Sigma) * ss) / df
    fact <- lgamma((df + p)/2) - lgamma(df/2) - 0.5 * (p * (log(pi) + log(df)) + sum(log(ev)))
    if (log)
        fact - 0.5 * (df + p) * log(1 + quad)
    else
        exp(fact) * ((1 + quad)^(- (df + p)/2))
}
