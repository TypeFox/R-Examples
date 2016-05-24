hot.csample<-function(X,g,method,n.simu,...)
 {

 g.levels<-levels(g)
 
 
 n<-dim(X)[1]
 p<-dim(X)[2]

 mean.total <- colMeans(X)
 
 Ts<-sweep(X,2,mean.total,"-")

 B<-crossprod(Ts)/n
 #B.inv<-solve(B)
 B.inv<-syminv(B)
 n.g<-by(g,g,length)
 
 T.g<-by(Ts,g,colMeans)
 
 Q.i<-sapply(T.g,my.quad.from,B.inv=B.inv,simplify=T)
 
 Q.2 <- sum(n.g*Q.i)
 names(Q.2) <- "Q.2"
 
 METHOD <- "Several samples location test using Hotellings T2"
 
 if (method=="approximation")
    { 
    parameter<-p*(nlevels(g)-1)
    names(parameter)<-"df"                       
    p.value<-1-pchisq(Q.2,parameter)
    }
    else
    {
    parameter <- n.simu
    names(parameter) <- "replications"
    statistics<-replicate(n.simu,perm.t(g[sample(1:n,n)],Ts,B.inv=B.inv, n.g=n.g ))
    p.value<- mean(statistics>Q.2)
    }
 names(Q.2) <- "Q.2"
 list(statistic=Q.2,p.value=p.value,parameter=parameter, method=METHOD)
 }


my.quad.from <- function (x,B.inv)
    {
    as.numeric(crossprod(x,B.inv)%*%x)
    }

perm.t <- function(g,Ts,B.inv,n.g)
    {
    T.g<-by(Ts,g,colMeans)
    Q.i<-sapply(T.g,my.quad.from,B.inv=B.inv,simplify=T)
    Q.2 <-sum(n.g*Q.i)
    Q.2
    }


CssTestOuter <- function(X,g,method,n.simu,...)
    {
    g.levels<-levels(g)
    n<-dim(X)[1]
    p<-dim(X)[2]
    
    Ts<-spatial.sign(X,center=TRUE, shape=FALSE,...)
    
    B<-crossprod(Ts)/n
    #B.inv<-solve(B)
    B.inv<-syminv(B)
    n.g<-by(g,g,length)
 
    T.g<-by(Ts,g,colMeans)
    
    Q.i<-sapply(T.g,my.quad.from,B.inv=B.inv,simplify=T)
 
    Q.2 <- sum(n.g*Q.i)
    names(Q.2) <- "Q.2"
    
    METHOD <- "Several samples location test using spatial signs"
 
 
    if (method=="approximation")
        { 
        parameter<-p*(nlevels(g)-1)
        names(parameter)<-"df"                       
        p.value<-1-pchisq(Q.2,parameter)
        }
    else
        {
        parameter <- n.simu
        names(parameter) <- "replications"
        statistics<-replicate(n.simu,perm.t(g[sample(1:n,n)],Ts,B.inv=B.inv, n.g=n.g ))
        p.value<-mean(statistics>Q.2)
        }
    list(statistic=Q.2,p.value=p.value,parameter=parameter, method=METHOD)
    }
    
    
 CsrTestOuter <- function(X,g,method,n.simu,...)
    {
    g.levels<-levels(g)
    n<-dim(X)[1]
    p<-dim(X)[2]
    
    Ts<-spatial.rank(X, shape=FALSE,...)
    
    B<-crossprod(Ts)/n
    #B.inv<-solve(B)
    B.inv<-syminv(B)
    n.g<-by(g,g,length)
 
    T.g<-by(Ts,g,colMeans)
    
    Q.i<-sapply(T.g,my.quad.from,B.inv=B.inv,simplify=T)
 
    Q.2 <- sum(n.g*Q.i)
    
    
    METHOD <- "Several samples location test using spatial ranks"
 
 
    if (method=="approximation")
        { 
        parameter<-p*(nlevels(g)-1)
        names(parameter)<-"df"                       
        p.value<-1-pchisq(Q.2,parameter)
        }
    else
        {
        parameter <- n.simu
        names(parameter) <- "replications"
        statistics<-replicate(n.simu,perm.t(g[sample(1:n,n)],Ts,B.inv=B.inv, n.g=n.g ))
        p.value<-mean(statistics>Q.2)
        }
    names(Q.2) <- "Q.2"
    list(statistic=Q.2,p.value=p.value,parameter=parameter, method=METHOD)
    }
    


CssTestInner <- function(X,g,method,n.simu,...)
    {
    g.levels<-levels(g)
    n<-dim(X)[1]
    p<-dim(X)[2]
    
    Ts<-spatial.sign(X,center=TRUE, shape=TRUE,...)

    n.g<-by(g,g,length)
 
    T.g<-by(Ts,g,colMeans)
    
    Q.i<-  sum(n.g*SpatialNP:::norm(matrix(unlist(T.g),ncol=p,byrow=T))^2)

 
    Q.2 <- p*Q.i
    
 
    METHOD <- "Equivariant several samples location test using spatial signs"
    
    if (method=="approximation")
        { 
        parameter<-p*(nlevels(g)-1)
        names(parameter)<-"df"                       
        p.value<-1-pchisq(Q.2,parameter)
        }
    else
        {
        parameter <- n.simu
        names(parameter) <- "replications"
        statistics<-replicate(n.simu,  sum(n.g*SpatialNP:::norm(matrix(unlist(by(Ts,g[sample(1:n,n)],colMeans)),ncol=p,byrow=T))^2))
        p.value<-mean(statistics>Q.i)
        }
    names(Q.2) <- "Q.2"
    list(statistic=Q.2,p.value=p.value,parameter=parameter, method=METHOD)
    }
    
    
 CsrTestInner <- function(X,g,method,n.simu,...)
    {
    g.levels<-levels(g)
    n<-dim(X)[1]
    p<-dim(X)[2]
    
    Ts<-spatial.rank(X, shape=TRUE,...)
  
    n.g<-by(g,g,length)
 
    T.g<-by(Ts,g,colMeans)
    
   
    Q.total<- mean(SpatialNP:::norm(Ts)^2)
    
    Q.i<-  sum(n.g*SpatialNP:::norm(matrix(unlist(T.g),ncol=p,byrow=T))^2)
 
    Q.2 <- p*Q.i / Q.total
    names(Q.2) <- "Q.2"
    
    METHOD <- "Equivariant several samples location test using spatial ranks"
 
    if (method=="approximation")
        { 
        parameter<-p*(nlevels(g)-1)
        names(parameter)<-"df"                       
        p.value<-1-pchisq(Q.2,parameter)
        }
    else
        {
        parameter <- n.simu
        names(parameter) <- "replications"
        statistics<-replicate(n.simu, sum(n.g*SpatialNP:::norm(matrix(unlist(by(Ts,g[sample(1:n,n)],colMeans)),ncol=p,byrow=T))^2))
        p.value<- mean(statistics>Q.i)
        }
    list(statistic=Q.2,p.value=p.value,parameter=parameter, method=METHOD)
    }
    
