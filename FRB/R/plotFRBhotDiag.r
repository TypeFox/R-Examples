diagplot.FRBhot <- function(x, ...) {

  FRBres <- x
  X <- FRBres$X
  n1 <- nrow(X)
  q <- ncol(X)

  if (is.null(x$Mu1)) { # one-sample test
  # compute the robust distances:
  mahD <- sqrt(mahalanobis(X, FRBres$Mu, FRBres$Sigma))

    cutoff <- sqrt(qchisq(.975, q))
    ylimm <- c(0,max(max(mahD), 2*cutoff))
    plot(1:n1, mahD, xlab="Index", ylab="Robust distance", cex.lab=1.3, pch=20, ylim=ylimm,cex.axis=1.1)
    title("Diagnostic plot based on robust estimates")
    abline(h=cutoff)

    # identify outliers
    inds <- 1:n1
    badinds <- inds[mahD > sqrt(qchisq(.975, q))]
    if (length(badinds) <= 20) {
      for (j in badinds)
        text(inds[j], mahD[j], j, pos=4, cex=0.8, offset=0.3)
    }
  }
  else { # two-sample test
    Y <- FRBres$Y
    n2 <- nrow(Y)
    n <- n1 + n2
      # compute the robust distances:
    mahDX <- sqrt(mahalanobis(X, FRBres$Mu1, FRBres$Sigma))
    mahDY <- sqrt(mahalanobis(Y, FRBres$Mu2, FRBres$Sigma))
    mahD <- c(mahDX, mahDY)
    
    gap <- n1/10
    indsplot <- c(1:n1, n1+n1/10+(1:n2))
    cutoff <- sqrt(qchisq(.975, q))
    ylimm <- c(0,max(max(mahD), 2*cutoff))
    plot(indsplot, mahD, xlab="sample 1 | sample 2", ylab="Robust distance", cex.lab=1.3, pch=20, ylim=ylimm, xaxt="n",cex.axis=1.1)
    title("Diagnostic plot based on robust estimates")
    axis(side=1, at = c(1,n1, n1+n1/10+1, n1+n1/10+n2), labels = c(1,n1,1,n2), cex.axis=1.1)
    abline(h=cutoff)
    abline(v = n1+(n1/10+1)/2, lwd=2, lty="dashed", col="red")

    # identify outliers
    indsX <- 1:n1
    badindsX <- indsX[mahDX > cutoff]
    if (length(badindsX) <= 12) {
      for (j in badindsX)
        text(indsX[j], mahDX[j], j, pos=4, cex=0.8, offset=0.3)
    }
    indsY <- 1:n2
    badindsY <- indsY[mahDY > cutoff]
    if (length(badindsY) <= 12) {
      for (j in badindsY)
        text(indsplot[n1+j], mahDY[j], j, pos=4, cex=0.8, offset=0.3)
    }

  }
    
}