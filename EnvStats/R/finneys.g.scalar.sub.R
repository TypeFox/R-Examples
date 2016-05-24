finneys.g.scalar.sub <-
function (m, z, n.terms) 
{
    if (n.terms <= 2) 
        stop("n.terms must be greater than 2")
    n <- n.terms - 1
    p <- 2:n
    num <- c(log(1), log(m) + log(abs(z)), 2 * p * log(m) + log(m + 
        2 * p) + p * log(abs(z)))
    denom <- c(log(1), log(m + 1), log(m) + log(m + 2) + cumsum(log(m + 
        2 * p)) + p * log(m + 1) + cumsum(log(p)))
    if (z > 0) 
        exp(num - denom)
    else sign((-1)^(0:n)) * exp(num - denom)
}
