### multinormal test based on kurtosis measures
### several options for p-value computation
### result of class htest

mvnorm.kur.test <- function(X, method = "integration", n.simu = 1000, na.action = na.fail)
    {
    DNAME<-deparse(substitute(X))
    
    X<-na.action(X)
    X<-as.matrix(X)
    
    n<-dim(X)[1]
    p<-dim(X)[2]
    
    method <- match.arg(method, c("integration", "satterthwaite", "simulation"), several.ok = FALSE)
  
    # internal function for test statistic computation
        
    .W.stat.func<-function(X)
        {
        C.1<-cov(X)
        C.2<-cov4(X)
        
        C.1.eigen <- eigen(C.1, symmetric = TRUE)
        C.1.sqrt.inv <- C.1.eigen$vectors %*% (diag(C.1.eigen$values^-0.5))%*%t(C.1.eigen$vectors)
        # better test statistic than in the paper since it guaranties affine invarince:
        frobenius.norm(C.1.sqrt.inv %*% C.2 %*% C.1.sqrt.inv - diag(p))^2
        }
    
    
    W.stat <- .W.stat.func(X)
    
    dfs <- c(0.5*p*(p+1)-1,1)
    chi.fac <- c(4*(p+4)/((p+2)^2),8/(p+2))
    p.value <- 1-pchisqsum(n*W.stat,chi.fac,dfs)
    
    # p-value computation based on the different options
    
    p.value<-switch(method, "integration"=1-pchisqsum(n*W.stat,df=dfs,a=chi.fac,method="integration"),
                  "satterthwaite"=1-pchisqsum(n*W.stat,df=dfs,a=chi.fac,method="satterthwaite"),
                  "simulation"={Ws <-replicate(n.simu,.W.stat.func(rmvnorm(n,rep(0,p))));
                  mean(Ws>W.stat)})
    
    
    
    STATISTIC<-n*W.stat
    names(STATISTIC)<-"W"
    METHOD<-"Multivariate Normality Test Based on Kurtosis"
    
    if(method=="simulation")
    {
    PARAMETER<-n.simu
    names(PARAMETER)<-c("replications")
    }
    else
    {
    PARAMETER<-c(chi.fac[1],dfs[1],chi.fac[2],dfs[2])
    names(PARAMETER)<-c("w1","df1","w2","df2")
    }
    
    PVAL<-p.value
    res<-list(method=METHOD,statistic=STATISTIC,data.name=DNAME,parameter=PARAMETER,p.value=PVAL)
    class(res)<-"htest"
    return(res)
    
    }
