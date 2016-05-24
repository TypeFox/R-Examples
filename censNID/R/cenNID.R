"cenNID" <-
  function(y, L=rep(-Inf, length(y)), U=rep(Inf,length(y)))
  {
    n <- length(y)
    x2 <- indP <- numeric(n)
    x1 <- y
    indP[L > -Inf] <- 1
    indP[U < Inf] <- -1
    iB <- (L > -Inf) & (U < Inf)
    if (any(iB)) {
        x1[iB]<-L[iB]
        x2[iB]<-U[iB]
        indP[iB] <- 2
    }
    xmean <- mean(y)
    xsigma <- sd(y)
    e1 <- e2 <- 10^(-6)
    cov <- matrix(numeric(4), ncol=2, nrow=2)
    maxits <- 100
    k <- 0
    nobs <- c(0,0,0,0)
    ifault <- 0
    outF <- .Fortran("em", as.integer(n), as.double(x1), as.double(x2), 
                     as.integer(indP), as.double(xmean), as.double(xsigma), 
                     as.double(e1), as.double(e2), as.integer(maxits),
                     as.double(cov), as.integer(nobs), as.integer(k), as.integer(ifault))
    muHat <- outF[[5]]
    sigHat <- outF[[6]]
    sds <- sqrt((outF[[10]])[c(1,4)])
    est <- matrix(c(muHat,sigHat,sds), ncol=2, nrow=2)
    dimnames(est) <- list(c("mean", "sd"), c("mle", "se(mle)"))
    covmat <- matrix(outF[[10]], ncol=2, nrow=2)
    dimnames(covmat) <- list(c("mean", "sd"), c("mean", "sd"))
    ans <- list(est=est, CovMat=covmat, nobs=outF[[11]], iterCount=outF[[12]], ifault=outF[[13]])
    ans
  }