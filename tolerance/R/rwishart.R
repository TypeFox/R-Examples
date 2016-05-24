#Internal Function

rwishart <- function (df, p) 
{
    X <- matrix(0, p, p)
    diag(X) <- sqrt(rchisq(p, df:(df - p + 1)))
    if (p > 1) {
        pseq <- 1:(p - 1)
        X[rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] = rnorm(p * 
            (p - 1)/2)
    }
    crossprod(X)
}
