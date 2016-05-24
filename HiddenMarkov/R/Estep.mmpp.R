"Estep.mmpp" <-
function (tau, Q, delta, lambda, fortran=TRUE) 
{
    #   scaled version of Ryden (1996)
    #   note that tau has already been differenced in BaumWelch.mmpp
    #   IT SHOULD BE CALLED SOMETHING DIFFERENT, SAY dtau
    m <- ncol(Q)
    n <- length(tau)
    x <- forwardback.mmpp(tau, Q, delta, lambda)
    logbeta <- x$logbeta
    logalpha <- x$logalpha
    alpha <- exp(logalpha)
    beta <- exp(logbeta)
    #   eigenvalue decomposition
    d <- x$eigenval
    S <- x$S
    Sinv <- x$Sinv
    diff <- outer(d, d, FUN="-")+diag(m)
    post0 <- Sinv %*% diag(lambda)
    A <- matrix(as.double(0), nrow=m, ncol=m)
    TT <- array(as.double(0), dim=c(n, m, m))
    if (fortran!=TRUE){
        #  loop5 using R code
        for (i in 1:n){
            expdtau <- exp(d*tau[i])
            difftau <- outer(expdtau, expdtau, FUN="-")
            TT[i,,] <- (difftau + diag(tau[i]*expdtau))/
                        diff/exp(x$scalefac[i])
        }
    } else{
        if (!is.double(d)) stop("d is not double precision")
        if (!is.double(tau)) stop("tau is not double precision")
        memory0 <- rep(as.double(0), m)
        loop5 <- .Fortran("loop5", m, n, d, tau, x$scalefac, diff, TT,
                          memory0, PACKAGE="HiddenMarkov")
        TT <- loop5[[7]]
    }
    if (fortran!=TRUE){
        #  loop6 using R code
        for (k in 1:m){
            pre <- S %*% diag(Sinv[,k])
            for (j in 1:m){
                post <- diag(S[j,]) %*% post0
                for (i in 1:n){
                    A[k,j] <- A[k,j] + alpha[i,] %*% pre %*% TT[i,,] %*%
                                       post %*% beta[(i+1),]
                }
            }
        }
    } else{
        memory0 <- matrix(as.double(0), ncol=m, nrow=m)
        memory1 <- matrix(as.double(0), ncol=m, nrow=m)
        memory2 <- matrix(as.double(0), ncol=m, nrow=m)
        memory3 <- rep(as.double(0), m)
        memory4 <- rep(as.double(0), m)
        loop6 <- .Fortran("loop6", m, n, TT, S, Sinv, post0, alpha,
                          beta, A, memory0, memory1, memory2, memory3,
                          memory4, PACKAGE="HiddenMarkov")
        A <- loop6[[9]]
    }
    if (n==1) B <- exp(logalpha[-1,]+logbeta[-1,])
    else B <- apply(exp(logalpha[-1,]+logbeta[-1,]), MARGIN=2, FUN=sum)
    return(list(A=A, B=B, logalpha=logalpha, logbeta=logbeta, LL=x$LL))
}

