### multinormal test based on skewness
### only assymptotic p-value returned
###

mvnorm.skew.test <- function(X, na.action = na.fail)
    {
    DNAME<-deparse(substitute(X))
    
    X<-na.action(X)
    X<-as.matrix(X)
    
    n<-dim(X)[1]
    p<-dim(X)[2]
    
    T.1<-colMeans(X)
    T.2<-mean3(X)
    C.M<-cov(X)
    
    U.stat <- crossprod((T.1-T.2), solve(C.M)) %*% (T.1-T.2)
    chi.fac <- 2*(p+2)/(p^2)
    p.value <- 1-pchisq(n*U.stat/chi.fac,p)
    
    STATISTIC<-n*U.stat/chi.fac
    names(STATISTIC)<-"U"
    METHOD<-"Multivariate Normality Test Based on Skewness"
    
    PARAMETER<-c(p)
    names(PARAMETER)<-c("df")
    PVAL<-p.value
    res<-list(method=METHOD, statistic=STATISTIC, data.name=DNAME, parameter=PARAMETER, p.value=PVAL)
    class(res)<-"htest"
    return(res)
    
    }
