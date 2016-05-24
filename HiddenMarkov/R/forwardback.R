"forwardback" <-
function(x, Pi, delta, distn, pm, pn = NULL, fortran=TRUE){
    m <- nrow(Pi)
    n <- length(x)
    dfunc <- makedensity(distn)
    prob <- matrix(as.double(0), nrow=n, ncol=m)
    for (k in 1:m)
        prob[,k] <- dfunc(x=x, getj(pm, k), pn, log=FALSE)
    #   forward probabilities alpha_ij
    phi <- as.double(delta)
    logalpha <- matrix(as.double(rep(0, m*n)), nrow=n)
    lscale <- as.double(0)
    if (fortran!=TRUE){
        #  loop1 using R code
        for (i in 1:n){
            if (i > 1) phi <- phi %*% Pi
            phi <- phi*prob[i,]
            sumphi <- sum(phi)
            phi <- phi/sumphi
            lscale <- lscale + log(sumphi)
            logalpha[i,] <- log(phi) + lscale
        }
        LL <- lscale
    } else{
        if (!is.double(Pi)) stop("Pi is not double precision")
        if (!is.double(prob)) stop("prob is not double precision")
        memory0 <- rep(as.double(0), m)
        loop1 <- .Fortran("loop1", m, n, phi, prob, Pi, logalpha,
                          lscale, memory0, PACKAGE="HiddenMarkov")
        logalpha <- loop1[[6]]
        LL <- loop1[[7]]
    }
    #   backward probabilities beta_ij
    logbeta <- matrix(as.double(rep(0, m*n)), nrow=n)
    phi <- as.double(rep(1/m, m))
    lscale <- as.double(log(m))
    if (fortran!=TRUE){
        #  loop2 using R code
        for (i in seq(n-1, 1, -1)){
            phi <- Pi %*% (prob[i+1,]*phi)
            logbeta[i,] <- log(phi) + lscale
            sumphi <- sum(phi)
            phi <- phi/sumphi
            lscale <- lscale + log(sumphi)
        }
    } else{
        memory0 <- rep(as.double(0), m)
        loop2 <- .Fortran("loop2", m, n, phi, prob, Pi, logbeta,
                          lscale, memory0, PACKAGE="HiddenMarkov")
        logbeta <- loop2[[6]]
    }
    return(list(logalpha=logalpha, logbeta=logbeta, LL=LL))
}

