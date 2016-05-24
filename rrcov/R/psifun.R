##  Internal functions, used for 
##  computing expectations
##

.chiInt <- function(p,a,c1)
##   partial expectation d in (0,c1) of d^a under chi-squared p
    return(exp(lgamma((p+a)/2)-lgamma(p/2))*2^{a/2}*pchisq(c1^2,p+a))

.chiInt2 <- function(p,a,c1)
##  partial expectation d in (c1,\infty) of d^a under chi-squared p
    return(exp(lgamma((p+a)/2)-lgamma(p/2))*2^{a/2}*(1-pchisq(c1^2,p+a)) )

##
##       derivatives of the above functions wrt c1
##
.chiIntD <- function(p,a,c1)
    return(exp(lgamma((p+a)/2)-lgamma(p/2))*2^{a/2}*dchisq(c1^2,p+a)*2*c1)

.chiInt2D <- function(p,a,c1)
    return(-exp(lgamma((p+a)/2)-lgamma(p/2))*2^{a/2}*dchisq(c1^2,p+a)*2*c1)


setMethod("iterM", "PsiFun", function(obj, x, t1, s, eps=1e-3, maxiter=20){
##  M-estimation from starting point (t1, s) for translated 
##  biweight with median scaling (Rocke & Woodruf (1993))
##  - obj contains the weight function parameters, e.g.
##      obj@c1 and obj@M are the constans for the translated 
##      biweight
##  - eps - precision for the convergence algorithm
##  - maxiter - maximal number of iterations
#
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    crit <- 100
    iter <- 1
    w1d <- w2d <- rep(1,n)
    while(crit > eps & iter <= maxiter) {
        t.old <- t1
        s.old <- s
        wt.old <- w1d
        v.old <- w2d

        ## compute the mahalanobis distances with the current estimates t1 and s
        d <- sqrt(mahalanobis(x, t1, s))
        
        # Compute k = sqrt(d[(n+p+1)/2])/M
        # and ajust the distances d = sqrt(d)/k 
        h <- (n+p+1)%/%2
        d <- d*sqrt(qchisq(h/(n+1), p))/(sort(d)[h])
        
        # compute the weights
        w1d <- wt(obj, d)
        w2d <- vt(obj, d)
        
        ## compute reweighted mean and covariance matrix
        t1 <- colSums(w1d*x)/sum(w1d)
        xx <- sqrt(w1d)*sweep(x, 2, t1)
        s <- p*(t(xx) %*% xx)/sum(w2d)

        ## check convergence
        crit <- max(abs(w1d-wt.old))/max(w1d)
        iter <- iter+1
    }
    return(list(t1=t1, s=s, iter=iter, wt=w1d, vt=w2d))
})
