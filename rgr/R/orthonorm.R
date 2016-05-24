orthonorm <-
function (p) 
{
    V <- matrix(0, nrow = p, ncol = p - 1)
    for (i in 1:(p - 1)) {
        V[1:i, i] <- 1/i
        V[i + 1, i] <- (-1)
        V[, i] <- V[, i] * sqrt(i/(i + 1))
    }
    return(V = V)
}
