monmlp.reshape <-
function(x, y, weights, hidden1, hidden2)
{
    N11 <- ncol(x)+1
    N12 <- hidden1
    N1 <- N11*N12
    W1 <- weights[1:N1]
    W1 <- matrix(W1, N11, N12)
    if (hidden2==0){
        N21 <- hidden1+1
        N22 <- ncol(y)
        N2 <- N1 + N21*N22
        W2 <- weights[(N1+1):N2]
        W2 <- matrix(W2, N21, N22)
        W.list <- list(W1=W1, W2=W2)
    } else{
        N21 <- hidden1+1
        N22 <- hidden2
        N2 <- N1 + N21*N22
        W2 <- weights[(N1+1):N2]
        W2 <- matrix(W2, N21, N22)
        N31 <- hidden2+1
        N32 <- ncol(y)
        N3 <- N2 + N31*N32
        W3 <- weights[(N2+1):N3]
        W3 <- matrix(W3, N31, N32)
        W.list <- list(W1=W1, W2=W2, W3=W3)
    }
    W.list
}

