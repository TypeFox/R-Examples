warfima <- function(z, order=c(0,0,0), lmodel=c("FD", "FGN", "PLA", "PLS", "NONE")) {
    lmodel <- match.arg(lmodel)
    p <- order[1]
    d <- order[2]
    q <- order[3]
    stopifnot(p >= 0, q >= 0, d >= 0)
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    stopifnot(is.wholenumber(p), is.wholenumber(d), is.wholenumber(q))
    w <- if(d>0) diff(z, differences=d) else z
    w <- w-mean(w)
    n <- length(w)
    alg <- 1
    binit <- numeric(p+q+1)
    binit[1] <- 1
    Ip <- (spec.pgram(w, fast=FALSE, detrend=FALSE, plot=FALSE, taper=0)$spec)/(2*pi)
    penaltyLoglikelihood <- (-n/2*log(sum(w^2)/n))*0.01
    Entropy<-function(beta, p, q) {
        alpha <- beta[1]
        phi <- theta <- numeric(0)
        if (alpha<=0 || alpha >=2 || (p>0||q>0)&&abs(beta[-1]) >= 1)
            EN <- -penaltyLoglikelihood else {
            if(p>0) phi   <- PacfToAR(beta[2:(p+1)])
            if(q>0) theta <- PacfToAR(beta[(p+2):(p+q+1)])
            fp <- sdfhd(n, alpha=alpha, phi=phi, theta=theta, lmodel=lmodel)
            sigHat <- mean(Ip/fp)
        EN <- ifelse(lmodel%in%c("FGN", "PLA", "PLS"), 2*sum(log(sigHat*fp)), 2*sigHat)
        }
        EN
    }
    if (p+q>0 || lmodel!="NONE") {
     ans<-optim(par=binit, fn=Entropy, p=p, q=q, method="L-BFGS-B", lower=c(0.01,rep(-0.99,p+q)), upper=c(1.99,rep(0.99,p+q)), control=list(trace=0))
     if(ans$convergence != 0) {#convergence problem. Use Nelder-Mead with penalty function
        alg<-2
        ans<-optim(par=binit, fn=Entropy,  p=p, q=q, method="Nelder-Mead")
        if(ans$convergence != 0) {#convergence problem. Use SANN with penalty function
            alg<-3
            ans<-optim(par=binit, fn=Entropy,  p=p, q=q, method="SANN")
            }
     }
     bHat <- ans$par
     wLL <- -ans$value
     convergence <- ans$convergence
    } else { #pure white noise
     bHat <- numeric(0)
     wLL <- -Entropy(1, 0, 0)
     convergence <- 0
    }
    alphaHat <- bHat[1]
    HHat <- 1-alphaHat/2
    dHat <- HHat - 0.5
    phiHat <- thetaHat <- numeric(0)
    if (p>0) phiHat <- PacfToAR(bHat[2:(p+1)])
    if (q>0) thetaHat <- PacfToAR(bHat[(2+p):(1+p+q)])
    list(bHat=bHat, alphaHat=alphaHat, HHat = HHat, dHat=dHat, phiHat=phiHat, thetaHat=thetaHat,
        LL=wLL, convergence=convergence, algorithm=alg)
}
