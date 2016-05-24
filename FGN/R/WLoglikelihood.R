WLoglikelihood <- function(z, beta, p, q, lmodel) {
    ialpha <- ifelse(identical(lmodel,"NONE"), 0, 1)
    alpha <- ifelse(ialpha==1, beta[1], 1)
    phi <- theta <- numeric(0)
    n <- length(z)
    Ip <- (spec.pgram(z, fast=FALSE, detrend=FALSE, plot=FALSE, taper=0)$spec)/(2*pi)
    if (alpha<=0 || alpha >=2 || (ialpha==1&&(p>0||q>0)&&abs(beta[-1]) >= 1))
        L <- LL <- (-n/2*log(sum(z^2)/n))-10^4 else {
          if(p>0) phi   <- PacfToAR(beta[(1+ialpha):(p+ialpha)])
          if(q>0) theta <- PacfToAR(beta[(ialpha+p+1):(p+q+ialpha)])
          fp <- sdfhd(n, alpha=alpha, phi=phi, theta=theta, lmodel=lmodel)
          sigHat <- mean(Ip/fp)
          L <- ifelse(lmodel%in%c("FGN", "PLA", "PLS"), -2*sum(log(sigHat*fp)), -2*sigHat)
          r <- switch(lmodel,
          FD=tacvfARFIMA(phi = phi, theta = theta, dfrac = 0.5-alpha/2, maxlag = n-1),
          FGN=tacvfARFIMA(phi = phi, theta = theta, H = 1-alpha/2, maxlag = n-1),
          PLA=tacvfARFIMA(phi = phi, theta = theta, alpha = alpha, maxlag = n-1),
          NONE=tacvfARFIMA(phi = phi, theta = theta, maxlag = n-1) )
          LL <- DLLoglikelihood(r, z)
     } 
     ans <- c(L, LL)
     names(ans) <- c("Whittle", "Exact")
     ans
}
