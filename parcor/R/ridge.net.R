`ridge.net` <-
function(X,lambda=NULL,plot.it=FALSE,scale=TRUE,k=10,verbose=FALSE){
    if (is.null(lambda)==TRUE){
        ss<-seq(-10,-1,length=1000)
        ss<-10^ss
        n<-nrow(X)
        nn<-n- floor(n/k)
        lambda<-ss*nn*ncol(X)
    }
    n <- nrow(X)
    p <- ncol(X)
    X <- scale(X,scale=scale)  # data needs to be centered and standardized
    B<- matrix(0, nrow=p, ncol=p) # matrix of regression coefficients
    lambda.opt<-rep(0,p)
    cat(paste("Performing local ridge regressions\n"))  
    cat(paste("Vertex no "))
  for (i in 1:p) ## visit all nodes
  { 
    if ((i/10)==floor(i/10)) {cat(paste(i,"..."))}
        noti <- (1:p)[-i]
        yi <- X[ , i]       # response
        Xi <- X[ , noti]    # predicted by all other nodes with i missing
        r.cv= ridge.cv(Xi,yi,lambda=lambda,scale=scale,plot.it=plot.it,k=k)
        B[i,-i]=r.cv$coefficients
        lambda.opt[i]=r.cv$lambda.opt
    }
    
    pcor<- Beta2parcor(B,verbose=verbose)  # compute matrix of partial correlations
  return(list(pcor=pcor))  
}
