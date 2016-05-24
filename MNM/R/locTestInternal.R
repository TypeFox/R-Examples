id.loc <- function(X, method,n.simu)
    {
    n<-dim(X)[1]
    p<-dim(X)[2]
    
    B<- crossprod(X)/n
    B.inv <- solve(B)
    Tstat<-colMeans(X)
    STATISTIC<-as.numeric( n * t(Tstat) %*% B.inv %*% Tstat)
    names(STATISTIC) <- "Q.2"
    
    f.stat<- function(Y,B.inv,n)  {Ts<-colMeans(Y)
                             as.numeric(n* t(Ts) %*% B.inv %*% Ts)
                            }
                            
    METHOD<-"Hotelling's one sample T2-test"
    
    res <- switch(method,
     "approximation" = {PVAL <- 1-pchisq(STATISTIC,df=p)
                        PARAMETER <- p
                        names(PARAMETER)<-"df"
                        list(statistic=STATISTIC,p.value=as.numeric(PVAL),method=METHOD,parameter=PARAMETER)}
     ,
     "signchange" = {statistics <- replicate (n.simu,f.stat(Y=sample(c(1,-1),n,replace=T) * X, B.inv=B.inv, n=n))
                    PVAL<-mean(statistics>STATISTIC)
                    PARAMETER <- n.simu
                    names(PARAMETER)<-"replications"
                    list(statistic=STATISTIC,p.value=as.numeric(PVAL),method=METHOD,parameter=PARAMETER)
                    })
    return(res)
    
    }

ssloc.outer<-function(X,method,n.simu)
    {
    # Assumes hypothesis is origin
    n<-dim(X)[1]
    p<-dim(X)[2]
    SCORES <- spatial.sign(X,center=FALSE,shape=FALSE)    
    B<- crossprod(SCORES)/n
    B.inv <- solve(B)
    Tstat<-colMeans(SCORES)
    STATISTIC<-as.numeric( n * t(Tstat) %*% B.inv %*% Tstat)
    names(STATISTIC) <- "Q.2"
    
    f.stat<- function(Y,B.inv,n)  {Ts<-colMeans(Y)
                             as.numeric(n* t(Ts) %*% B.inv %*% Ts)
                            }
    
    METHOD<-"One sample spatial sign test\n using outer standardization"
    
    
    res <- switch(method,
     "approximation" = {PVAL <- 1-pchisq(STATISTIC,df=p)
                        PARAMETER <- p
                        names(PARAMETER)<-"df"
                        list(statistic=STATISTIC,p.value=as.numeric(PVAL),method=METHOD,parameter=PARAMETER)}
     ,
     "signchange" = {statistics <- replicate (n.simu,f.stat(Y=sample(c(1,-1),n,replace=T) * SCORES, B.inv=B.inv, n=n))
                    PVAL<-mean(statistics>STATISTIC)
                    PARAMETER <- n.simu
                    names(PARAMETER)<-"replications"
                    list(statistic=STATISTIC,p.value=as.numeric(PVAL),method=METHOD,parameter=PARAMETER)
                    })
    return(res)
    
    }


ssloc.inner<-function(X,method,n.simu)
    {
    # Assumes hypothesis is origin
    n<-dim(X)[1]
    p<-dim(X)[2]
    SCORES <- spatial.sign(X,center=FALSE,shape=TRUE)    
    Tstat<- sum(colMeans(SCORES)^2)
    STATISTIC<- n * p * Tstat
    METHOD<-"One sample spatial sign test\n using inner standardization"
    names(STATISTIC) <- "Q.2"
    
    res <- switch(method,
     "approximation" = {PVAL <- 1-pchisq(STATISTIC,df=p)
                        PARAMETER <- p
                        names(PARAMETER)<-"df"
                        list(statistic=STATISTIC,p.value=as.numeric(PVAL),method=METHOD,parameter=PARAMETER)}
     ,
     "signchange" = {statistics <- replicate (n.simu,sum(colMeans(sample(c(1,-1),n,replace=T) * SCORES)^2))
                    PVAL <- mean(statistics>Tstat)
                    PARAMETER <- n.simu
                    names(PARAMETER)<-"replications"
                    list(statistic=STATISTIC,p.value=as.numeric(PVAL),method=METHOD,parameter=PARAMETER)
                    })
    
    return(res)
    }



srloc.outer<-function(X,method,n.simu)
    {
    # Assumes hypothesis is origin
    n<-dim(X)[1]
    p<-dim(X)[2]
    SCORES <- SpatialNP:::signranks(X)   
    B<- crossprod(SCORES)/n
    B.inv <- solve(B)
    Tstat<-colMeans(SCORES)
    STATISTIC<- as.numeric(n * t(Tstat) %*% B.inv %*% Tstat)
    METHOD<-"One sample spatial signed-rank test\n using outer standardization"
    names(STATISTIC) <- "Q.2"
    
    f.stat<- function(Y,B.inv,n)  {Ts<-colMeans(Y)
                             as.numeric(n* t(Ts) %*% B.inv %*% Ts)
                            }
    
    res <- switch(method,
     "approximation" = {PVAL <- 1-pchisq(STATISTIC,df=p)
                        PARAMETER <- p
                        names(PARAMETER)<-"df"
                        list(statistic=STATISTIC,p.value=as.numeric(PVAL),method=METHOD,parameter=PARAMETER)}
     ,
     "signchange" = {statistics <- replicate (n.simu,f.stat(Y=sample(c(1,-1),n,replace=T) * SCORES, B.inv=B.inv, n=n))
                    PVAL <- mean(statistics>STATISTIC)
                    PARAMETER <- n.simu
                    names(PARAMETER)<-"replications"
                    list(statistic=STATISTIC,p.value=as.numeric(PVAL),method=METHOD,parameter=PARAMETER)
                    }
     )
    return(res)
    
    }


srloc.inner<-function(X,method,n.simu)
    {
    # Assumes hypothesis is origin
    n<-dim(X)[1]
    p<-dim(X)[2]
    SCORES <- spatial.signrank(X,center=F)   
    Q1<- sum(colMeans(SCORES)^2)
    Q2<- mean(SpatialNP:::norm(SCORES)^2)
    ratio<-Q1 / Q2
    STATISTIC<- n * p * ratio
    METHOD<-"One sample spatial signed-rank test\n using inner standardization"
    names(STATISTIC) <- "Q.2"
    
     f.stat<- function(Y)  {Q1<- sum(colMeans(Y)^2)
                            Q2<- mean(SpatialNP:::norm(Y)^2)
                            ratio<-Q1 / Q2
                            }
    
    res <- switch(method,
     "approximation" = {PVAL <- 1-pchisq(STATISTIC,df=p)
                        PARAMETER <- p
                        names(PARAMETER)<-"df"
                        list(statistic=STATISTIC,p.value=as.numeric(PVAL),method=METHOD,parameter=PARAMETER)}
     ,
     "signchange" = {statistics <- replicate (n.simu,f.stat(sample(c(1,-1),n,replace=T) * SCORES))
                    PVAL <- mean(statistics>ratio)
                    PARAMETER <- n.simu
                    names(PARAMETER)<-"replications"
                    list(statistic=STATISTIC,p.value=as.numeric(PVAL),method=METHOD,parameter=PARAMETER)
                    }
    )
    return(res)
    }
