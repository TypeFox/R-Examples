# Solves 1/2n |Y-Xbeta|^2_F+lambda|beta|_1
.solveMyLassoF_c<-function(X,Y,lambda,binit=NULL,eps=1e-4,maxiter=1000){
    p=ncol(X)
    n=nrow(Y)
    r=ncol(Y)
    if (is.null(binit)){binit=matrix(0,p,r)}
    beta=binit
    niter=0
    
    out=.C("solveMyLassoF",as.double(X),as.double(Y),as.double(beta),as.double(lambda),as.integer(p),as.integer(n),as.integer(r),as.double(eps),as.integer(maxiter),as.integer(niter))
    
    if (out[[10]]==maxiter){warning("Failed to converge, try increasing the number of iterations")}
    return(matrix(out[[3]],p,r))
}