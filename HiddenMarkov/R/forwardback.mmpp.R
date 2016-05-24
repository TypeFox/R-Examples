"forwardback.mmpp" <-
function(tau, Q, delta, lambda, fortran=TRUE, fwd.only=FALSE){
    #   note that tau has already been differenced in BaumWelch.mmpp
    #   IT SHOULD BE CALLED SOMETHING DIFFERENT, SAY dtau
    m <- nrow(Q)
    n <- length(tau)
    Lambda <- diag(lambda)
    ##    eigenvalue decomposition
    decomp <- eigen(Q-Lambda, symmetric=FALSE)
    if (any(duplicated(decomp$values))) stop("repeated eigenvalues")
    S <- decomp$vectors
    Sinv <- solve(S)
    eigenval <- decomp$values
    ##  scaled forward probabilities
    phi <- as.double(delta)
    logalpha <- matrix(as.double(rep(0, m*(n+1))), nrow=(n+1))
    logalpha[1,] <- log(phi)
    scalefac <- as.double(rep(0, n))
    post <- Sinv %*% Lambda
    psi <- array(as.double(0), dim=c(n,m,m))
    if (fortran!=TRUE){
        #   loop3 using R code
        for (i in 1:n){
            psi[i,,] <- S %*% diag(exp(eigenval*tau[i])) %*% post
            phi <- phi %*% psi[i,,]
            sumphi <- sum(phi)
            scalefac[i] <- log(sumphi)
            phi <- phi/sumphi
            logalpha[i+1,] <- log(phi)
        }
    } else{
        memory0 <- matrix(as.double(0), ncol=m, nrow=m)
        memory1 <- matrix(as.double(0), ncol=m, nrow=m)
        memory2 <- rep(as.double(0), m)
        if (!is.double(S)) stop("Eigenvectors are not double precision")
        loop3 <- .Fortran("loop3", m, n, phi, S, eigenval, logalpha,
                          scalefac, tau, post, psi, memory0, memory1,
                          memory2,
                          PACKAGE="HiddenMarkov", NAOK=TRUE)
        logalpha <- loop3[[6]]
        psi <- loop3[[10]]
        scalefac <- loop3[[7]]
    }
    if (fwd.only)
        return(list(logalpha=logalpha, eigenval=eigenval, S=S, 
                    Sinv=Sinv, scalefac=scalefac, LL=sum(scalefac)))
    ##  scaled backward probabilities
    logbeta <- matrix(rep(as.double(0), m*(n+1)), nrow=(n+1))
    phi <- rep(as.double(1/m), m)
    if (fortran!=TRUE){
        #   loop4 using R code
        lscale <- log(m)
        logck <- 0
        for (i in seq(n, 1, -1)){
            phi <- psi[i,,] %*% phi
            logck <- logck + scalefac[i] 
            logbeta[i,] <- log(phi) + lscale - logck
            sumphi <- sum(phi)
            phi <- phi/sumphi
            lscale <- lscale + log(sumphi)
        }
    }else{
        memory0 <- matrix(as.double(0), ncol=m, nrow=m)
        memory1 <- rep(as.double(0), m)
        loop4 <- .Fortran("loop4", m, n, phi, logbeta, scalefac, psi,
                          memory0, memory1, PACKAGE="HiddenMarkov")
        logbeta <- loop4[[4]]
    }
    return(list(logalpha=logalpha, logbeta=logbeta,
                eigenval=eigenval, S=S, Sinv=Sinv,
                scalefac=scalefac, LL=sum(scalefac)))
}

