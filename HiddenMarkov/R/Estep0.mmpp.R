"Estep0.mmpp" <-
function (tau, Q, delta, lambda) 
{
    m <- ncol(Q)
    n <- length(tau)
    logbeta <- backward0.mmpp(tau, Q, lambda)
    logalpha <- forward0.mmpp(tau, Q, delta, lambda)
    #   eigenvalue decomposition
    decomp <- eigen(Q-diag(lambda), symmetric=FALSE)
    d <- decomp$values
    S <- decomp$vectors
    Sinv <- solve(S)
    diff <- outer(d, d, FUN="-")
    post0 <- Sinv %*% diag(lambda)
    A <- matrix(0, nrow=m, ncol=m)
    for (i in 1:m){
        pre <- S %*% diag(Sinv[,i])
        for (j in 1:m){
            post <- diag(S[j,]) %*% post0
            for (k in 1:n){
                difftau <- outer(exp(d*tau[k]), exp(d*tau[k]), FUN="-")
                TT <- (difftau + diag(tau[k]*exp(d*tau[k])))/
                            (diff+diag(m))
                A[i,j] <- A[i,j] + exp(logalpha[k,]) %*% pre %*% TT %*%
                                   post %*% exp(logbeta[(k+1),])
            }
        }
    }
    if (n==1) B <- exp(logalpha[-1,]+logbeta[-1,])
    else B <- apply(exp(logalpha[-1,]+logbeta[-1,]), MARGIN=2, FUN=sum)
    return(list(A=A, B=B, logalpha=logalpha, logbeta=logbeta))
}

