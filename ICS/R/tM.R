### Internal function for joint norm of a location and matrix
### used in tM to decide about convergence
###

.norm.mu.V<-function(a,B,A)
    {
    A.inv<-solve(A)
    BA.inv <- B %*% A.inv
    square <- t(a) %*% A.inv %*% (a) + sum(diag( BA.inv %*% BA.inv))
    as.vector(sqrt(square))
    }

### Internal function for tM
### computes location and scatter according to Algorithm 1
### in Kent and Tyler

.alg1<-function(X,mu.init,V.init,nu,eps,maxiter)
    {
    p<-dim(X)[2]
    n<-dim(X)[1]
    V.i<-V.init
    mu.i<-mu.init
    iter<-0
    differ <- Inf
    while(differ > eps)
        {
        iter<- iter+1
        u<- (nu+p)/(nu+ mahalanobis(X,mu.i,V.i))
        mu.new<-colMeans(sweep(X,1,u,"*"))/mean(u)
        X.center2 <- sweep(X,2,mu.new)
        V.new <- t(sweep(X.center2,1,u,"*")) %*% X.center2 /n
        differ<- .norm.mu.V(a=mu.new-mu.i, B=V.new-V.i, A=V.new)
        V.i<-V.new
        mu.i<-mu.new
        if (iter>= maxiter) stop("maxiter reached without convergence")
        }
    return(list(mu=mu.new, V=V.new, iter=iter))
    }

### Internal function for tM
### computes location and scatter according to Algorithm 2
### in Kent and Tyler

.alg2<-function(X,mu.init,V.init, gamma.init,nu,eps,maxiter)
    {
    p<-dim(X)[2]
    n<-dim(X)[1]
    V.i<-V.init
    mu.i<-mu.init
    gamma.i<-gamma.init
    iter<-0
    differ <- Inf
    while(differ > eps)
        {
        iter<- iter+1
        w<- (nu+p)/(nu -1 + 1/gamma.i + (1/gamma.i)*mahalanobis(X,mu.i,V.i))
        gamma.new <-mean(w)
        mu.new<-colMeans(sweep(X,1,w,"*"))/gamma.new
        X.center2 <- sweep(X,2,mu.new)
        V.new <- (t(sweep(X.center2,1,w,"*")) %*% X.center2 /n) /gamma.new
        differ<- .norm.mu.V(a=(mu.new-mu.i), B=(V.new-V.i), A=V.new)
        gamma.i<-gamma.new
        V.i<-V.new
        mu.i<-mu.new
        if (iter>= maxiter) stop("maxiter reached without convergence")
        }
    return(list(mu=mu.new, V=V.new,gam=gamma.i, iter=iter))
    }

### Internal function for tM
### computes location and scatter according to Algorithm 3
### in Kent and Tyler

.alg3<-function(X,mu.init,V.init, nu,eps,maxiter)
    {
    p<-dim(X)[2]
    n<-dim(X)[1]
    V.i<-V.init
    mu.i<-mu.init
    gamma.i<-1
    iter<-0
    differ <- Inf
    while(differ > eps)
        {
        iter<- iter+1
        w<- (nu+p)/(nu -1 + 1/gamma.i + (1/gamma.i)*mahalanobis(X,mu.i,V.i))
        gamma.new <-mean(w)
        mu.new<-colMeans(sweep(X,1,w,"*"))/gamma.new
        X.center2 <- sweep(X,2,mu.new)
        V.new <- (t(sweep(X.center2,1,w,"*")) %*% X.center2 /n) /gamma.new
        differ<- .norm.mu.V(a=mu.new-mu.i, B=V.new-V.i, A=V.new)
               V.i<-V.new
        mu.i<-mu.new
        if (iter>= maxiter) stop("maxiter reached without convergence")
        }
    return(list(mu=mu.new, V=V.new, iter=iter))
    }

### function to compute the location and scatter for a multivariate t-distribution
### for a given degree of freedom
### is basically a wrapper around the subroutines .alg1, .alg2 and .alg3

tM<-function(X,df=1,alg="alg3",mu.init=NULL,V.init=NULL,gamma.init=NULL,eps=1e-06,maxiter=100, na.action=na.fail)
    {
    X<-na.action(X)
    X<-as.matrix(X)
    if (is.null(mu.init)) mu.init<-colMeans(X)
    if (is.null(V.init)) V.init<-cov(X)
    
    alg <- match.arg(alg,c("alg1","alg2","alg3"))
    
    if (alg!="alg2") if (!is.null(gamma.init)) warning("A initial value for gamma is only for alg2 needed")
    if (is.null(gamma.init)) gamma.init<-1
    res <-  switch(alg, 
                "alg1"=.alg1(X,mu.init=mu.init,V.init=V.init,nu=df,eps=eps,maxiter=maxiter),
                "alg2"=.alg2(X,mu.init=mu.init,V.init=V.init,nu=df,gamma.init=gamma.init,eps=eps,maxiter=maxiter),
                "alg3"=.alg3(X,mu.init=mu.init,V.init=V.init,nu=df,eps=eps,maxiter=maxiter)
                )
    return(res)
    }
    
