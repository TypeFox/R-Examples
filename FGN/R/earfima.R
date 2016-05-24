earfima <- function(z, order=c(0,0,0), lmodel=c("FD", "FGN", "PLA", "NONE")) {
    lmodel <- match.arg(lmodel)
    p <- order[1]
    d <- order[2]
    q <- order[3]
    stopifnot(p >= 0, q >= 0, d >= 0)
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)
          abs(x - round(x)) < tol
    stopifnot(is.wholenumber(p), is.wholenumber(d), is.wholenumber(q))
    w <- if(d>0) diff(z, differences=d) else z
    w <- w-mean(w)
    n <- length(w)
    alg <- 1
    binit <- numeric(p+q+1)
    binit[1] <- 1
    penaltyLoglikelihood <- (-n/2*log(sum(w^2)/n))*0.01
    Entropy<-function(beta, p, q) {
        alpha <- beta[1]
        phi <- theta <- numeric(0)
        if (alpha<=0 || alpha >=2 || (p>0||q>0)&&abs(beta[-1]) >= 1)
            LL <- -penaltyLoglikelihood else {
          if(p>0) phi   <- PacfToAR(beta[2:(p+1)])
          if(q>0) theta <- PacfToAR(beta[(p+2):(p+q+1)])
          r <- switch(lmodel,
            FD=tacvfARFIMA(phi = phi, theta = theta, dfrac = 0.5-alpha/2, maxlag = n-1),
            FGN=tacvfARFIMA(phi = phi, theta = theta, H = 1-alpha/2, maxlag = n-1),
            PLA=tacvfARFIMA(phi = phi, theta = theta, alpha = alpha, maxlag = n-1),
            NONE=tacvfARFIMA(phi = phi, theta = theta, maxlag = n-1) )
#needed in some cases, eg. NileMin with p=1, q=3, lmodel="FGN"
          LL <- try(-DLLoglikelihood(r, w), silent=TRUE)
          LL <- ifelse(is.numeric(LL), LL, -penaltyLoglikelihood)
          }
        LL
        }
    if (p+q>0 || lmodel!="NONE") {
     ans<-optim(par=binit, fn=Entropy, p=p, q=q, method="L-BFGS-B",
        lower=c(0.01,rep(-0.99,p+q)), upper=c(1.99,rep(0.99,p+q)), control=list(trace=0))
     if(ans$convergence != 0) {#convergence problem. Use Nelder-Mead with penalty function
        alg<-2
        ans<-optim(par=binit, fn=Entropy,  p=p, q=q, method="Nelder-Mead")
        if(ans$convergence != 0) {#convergence problem. Use SANN with penalty function
            alg<-3
            ans<-optim(par=binit, fn=Entropy,  p=p, q=q, method="SANN")
            }
     }
     bHat <- ans$par
     LL <- -ans$value
     convergence <- ans$convergence
    } else {
     bHat <- numeric(0)
     LL <- -Entropy(1, 0, 0)
     convergence <- 0
    }
    alphaHat <- bHat[1]
    HHat <- 1-alphaHat/2
    dHat <- HHat - 0.5
    phiHat <- thetaHat <- numeric(0)
    list(bHat=bHat, alphaHat=alphaHat, HHat = HHat, dHat=dHat, phiHat=phiHat, thetaHat=thetaHat,
       LL=LL, convergence=convergence, algorithm=alg)
}
