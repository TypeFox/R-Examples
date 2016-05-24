cadence.reshape <-
function(x, weights, n.hidden, distribution)
{
    N11 <- ncol(x)+1
    N12 <- n.hidden
    N1 <- N11*N12
    W1 <- weights[1:N1]
    W1 <- matrix(W1, N11, N12)
    N21 <- n.hidden+1
    N22 <- length(distribution$parameters)
    N2 <- N1 + N21*N22
    W2 <- weights[(N1+1):N2]
    W2 <- matrix(W2, N21, N22)
    list(W1=W1, W2=W2)
}

