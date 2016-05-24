# Solves 1/2n |Y-Xbeta|^2_2+lambda|beta|_1
.solveMyLasso_c<-function(X,Y,lambda,binit=NULL,eps=1e-4,maxiter=1000){
    p=ncol(X)
    n=length(Y)
    if (is.null(binit)){binit=rep(0,p)}
    beta=binit
    niter=0
    
    out=.C("solveMyLasso",as.double(X),as.double(Y),as.double(beta),as.double(lambda),as.integer(p),as.integer(n),as.double(eps),as.integer(maxiter),as.integer(niter))
    
    if (out[[9]]==maxiter){warning("Failed to converge, try increasing the number of iterations")}
    return(out[[3]])
}