.logRhoCdf <- function(y, lambda, m, A, k, n)
{
    index <- (y <= A)
    A.i <- A[index]
    k.i <- k[index]
    term.j <- 0
    for (j in 0:m)
        term.j <- term.j + stats::pgamma(q = log(A.i / y) * lambda / k.i,
                                         shape = j + 1,
                                         rate = 1,
                                         lower.tail = FALSE)

    term.i <- term.j / (m + 1)
    return(log((sum(term.i) / n)))
}
