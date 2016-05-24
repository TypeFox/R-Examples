`pls.net` <-
function(X,scale=TRUE,k=10,ncomp=15,verbose=FALSE){
    p=ncol(X)
    n=nrow(X)
    if (is.null(ncomp)){
        ncomp = min(n - 1, ncol(X))
    }
    k=floor(k)
    k=max(1,k)
    if (k>n){
        cat(paste("k exceeds the number of observations. Leave-one-out is applied.\n"))
        k=n
    }
    B=matrix(0,p,p) # regression coefficients
    m=vector(length=p) # optimal number of components
    cat(paste("Performing local pls regressions\n"))
    kernel=FALSE
    if (n<(p-1)){
        kernel=TRUE
    }
    cat(paste("Vertex no "))
  for (i in 1:p) ## visit all nodes
  { 
    if ((i/10)==floor(i/10)) {cat(paste(i,"..."))}
        Xi=X[,-i]
        yi=X[,i]
        fit=penalized.pls.cv(Xi,yi,scale=scale,k=k,ncomp=ncomp)
        #fit<-penalized.pls.cv(Xi,yi,P=NULL,scale=scale,k=k,ncomp=ncomp,kernel=kernel)
        B[i,-i]=fit$coefficients
        m[i]=fit$ncomp.opt
    }
    cat(paste("\n"))
    pcor <- Beta2parcor(B,verbose=verbose)
    return(list(pcor=pcor,m=m))
}
