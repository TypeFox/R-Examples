MardiaTest <-
function(data)
  {
    #SU: 06/15/2006: the following works in R only if you have the fSeries package loaded.
    # It loads the RMetrics fCalendar timeSeries class.  Since the other fixes for timeSeries
    # assume we are using the RMetrics package, this is OK.

    d <- dim(data)[2]
    n <- dim(data)[1]
    Xbar <- apply(data,2,mean)
    Xbar.matrix <- matrix(Xbar,nrow=n,ncol=d,byrow=TRUE)
    standardised <- data -Xbar.matrix
    S <- var(data)
    A <- t(chol(S))
    Ainv <- solve(A)
    Zdata <- Ainv %*% t(standardised)
    Zdata <- t(Zdata)
    D2.check <- mahalanobis(data,Xbar,S)
    Dij <- Zdata %*% t(Zdata)
    D2 <- diag(Dij)
    K3 <- mean(as.vector(Dij)^3)
    K4 <- mean(D2^2)
    statK3 <- n*K3/6
    df <- d*(d+1)*(d+2)/6
    K3.pval <- 1-pchisq(statK3,df)
    mn <- d*(d+2)
    vr <- 8*d*(d+2)/n
    K4.stat <- (K4-mn)/sqrt(vr)
    K4.pval <- 1-pnorm(abs(K4.stat))
    c(K3,K3.pval,K4,K4.pval)

outList = list ("skewness" = K3, "p.value" = K3.pval, "kurtosis" = K4, "p.value" = K4.pval)
 return(outList)
	}
