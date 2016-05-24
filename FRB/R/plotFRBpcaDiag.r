
diagplot.FRBpca <- function(x, EIF=TRUE, ...) {

  FRBres <- x
  Y <- FRBres$Y
  n <- nrow(Y)
  q <- ncol(Y)
  # compute the robust distances:
  mahD <- sqrt(mahalanobis(Y, FRBres$est$Mu, FRBres$est$Sigma))
  
  if (EIF) {
    detSigma <- det(FRBres$est$Sigma)^(-1/q)

    IFnorms <- rep(0,n)
    for (obs in 1:n) {

      IFGamma <- detSigma * mahD[obs]^2 * ( crossprod(Y[obs,]-FRBres$est$Mu)/(mahD[obs]^2) - FRBres$est$Sigma/q )
      IFeigvecs <- matrix(0,q,q)
      for (i in 1:q) {
          IFeigveci <- rep(0,q)
          for (j in 1:q) {
              if (j!=i) IFeigveci <- IFeigveci + 1/(FRBres$eigval[i]-FRBres$eigval[j])*(t(FRBres$eigvec[,j])%*%IFGamma%*%FRBres$eigvec[,i]) * FRBres$eigvec[,j]
          }
          IFeigvecs[,i] <- IFeigveci
      }
      IFnorms[obs] <- sqrt( 1/q^2 * sum(IFeigvecs^2) )
    }

    # cutoff by simulation
    nsim <- 50
    simIFnorms <- matrix(0,nsim,n);

    # sqrt for Sigma matrix:
    eig <- eigen(FRBres$est$Sigma)
    V <- eig$vectors
    sqrtSigma <- V %*% diag(sqrt(eig$values)) %*% t(V)
    simMu <- FRBres$est$Mu

    for (iter in 1:nsim) {
        simY <- matrix(rnorm(n*q),ncol=q) %*% sqrtSigma + matrix(rep(simMu,n),n,byrow=TRUE)

        clasSigma <- cov(simY)
        clasGamma <- det(clasSigma)^(-1/q) * clasSigma
        clasMu <- apply(simY, 2, mean)
        simmahD <- sqrt(mahalanobis(simY, clasMu, clasSigma))

        # compute corresponding eigenvalue/eigenvector estimates
        eigres <- eigen(clasGamma)
        IXGamma <- order(eigres$values, decreasing=TRUE)
        eigvals <- eigres$values[IXGamma]
        eigvecs <- eigres$vectors[,IXGamma]
        detSigma <- det(clasSigma)^(-1/q)

        for (obs in 1:n) {
          IFGamma <- detSigma * simmahD[obs]^2 * ( tcrossprod(simY[obs,]-clasMu)/(simmahD[obs]^2) - clasSigma/q )
          IFeigvecs <- matrix(0,q,q)
          for (i in 1:q) {
              IFeigveci <- rep(0,q)
              for (j in (1:q)[-i])
                IFeigveci <- IFeigveci + 1/(eigvals[i]-eigvals[j])*(t(eigvecs[,j])%*%IFGamma%*%eigvecs[,i]) * eigvecs[,j]
              IFeigvecs[,i] <- IFeigveci
          }
          simIFnorms[iter,obs] <- sqrt( 1/q^2 * sum(IFeigvecs^2) )
        }
    }

    cutoff <- quantile(unlist(simIFnorms),0.95)

    xlimm <- range(IFnorms)
    xlimm[2] <- xlimm[2] + (xlimm[2]-xlimm[1])*0.03
    plot(IFnorms, mahD, xlab="Empirical influence", ylab="Robust distance", cex.lab=1.3, pch=20,xlim=xlimm)
    title("Diagnostic plot based on robust estimates")
    abline(h=sqrt(qchisq(.975, q)))
    abline(v=cutoff)

    # identify outliers
    inds <- 1:n
    badinds <- inds[mahD > sqrt(qchisq(.975, q))]
    if (length(badinds) <= 20) {
      for (j in badinds)
        text(IFnorms[j], mahD[j], j, pos=4, cex=0.8, offset=0.3)
    }
  }
  else {
    plot(1:n, mahD, xlab="Index", ylab="Robust distance", cex.lab=1.3, pch=20)
    title("Diagnostic plot based on robust estimates")
    abline(h=sqrt(qchisq(.975, q)))

    # identify outliers
    inds <- 1:n
    badinds <- inds[mahD > sqrt(qchisq(.975, q))]
    if (length(badinds) <= 20) {
      for (j in badinds)
        text(inds[j], mahD[j], j, pos=4, cex=0.8, offset=0.3)
    }

  }
}
