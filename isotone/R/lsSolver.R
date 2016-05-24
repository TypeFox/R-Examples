#Least squares with diagonal weights
#returns fitted values y, Lagrange multiplier (lambda) lbd, function value f, gradient gy

lsSolver <- function(z, a, extra) 
{
    x <- z
    if ((is.null(extra$weights)) || (is.null(extra$y))) stop("lsSolver needs the additional arguments y and weights!")
    w <- extra$weights                           #weights
    z <- extra$y                                 #response
    n <- length(z)
    if (length(a)==0) return(list(x=z,lbd=NULL,f=0,gx=rep(0,n))) #no active set, break
    if (is.vector(a)) a <- matrix(a,1,length(a)) #only 1 active set (as matrix)
    
    indi <- mkIndi(a,n)                          #compute indicators
    h <- crossprod(indi, w*indi)
    r <- drop(crossprod(indi,w*z))
    b <- solve(h,r)
    y <- drop(indi%*%b)
    gy <- 2*w*(y-z)                              #gradient
    lbd <- mkLagrange(a, gy)
    f <- sum(w*(y-z)^2)                          #value target function
    return(list(x = y, lbd = lbd, f = f, gx = gy))
}
