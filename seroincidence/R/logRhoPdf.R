.logRhoPdf <- function(y, lambda, m, A, k, n)
{
    index <- (y <= A)
    A.i <- A[index]
    k.i <- k[index]
    delta_t.p <- (m + 1)/lambda
    term.i <- stats::pgamma(q = log(A.i / y)/ k.i,
                            shape = m + 1,
                            rate = lambda,
                            lower.tail = FALSE) / (k.i * delta_t.p)

    return(log((sum(term.i) / n) / y))
}
