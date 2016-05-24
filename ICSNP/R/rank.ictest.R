`rank.ictest`<-function(X,...)
    {
    UseMethod("rank.ictest")
    }

`rank.ictest.default` <-
function(X,mu=NULL,scores="rank",method = "approximation", n.simu = 1000,na.action=na.fail,...)
    {
    DNAME<-deparse(substitute(X))
    
    p<-dim(X)[2]

    if (is.null(mu)) mu<-rep(0,p)
    else if (length(mu)!=p) stop("length of 'mu' must equal the number of columns of 'X'")
    
    X<-na.action(X)
    n<-dim(X)[1]
    
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric")
    X<-as.matrix(X)

    X<-sweep(X,2,mu)
    score<-match.arg(scores,c("sign","rank","normal"))
    method <- match.arg(method,c("approximation","simulation","permutation"))
    

    res1<-switch(score,
        "sign"={
               test<-Q.Test(X,score)
               STATISTIC<-test$test.statistic
               names(STATISTIC)<-"Q.S"
               PVAL<-test$p.value
               METHOD<-"MARGINAL SIGN TEST (assuming the independent component model)"
               RVAL<-list(statistic=STATISTIC,p.value=PVAL,method=METHOD)
               
               RVAL}
        ,
        "rank"={
               test<-Q.Test(X,score)
               STATISTIC<-test$test.statistic
               names(STATISTIC)<-"Q.W"
               PVAL<-test$p.value
               METHOD<-"MARGINAL SIGNED RANK TEST (assuming the independent component model)"
               RVAL<-list(statistic=STATISTIC,p.value=PVAL,method=METHOD)
               
               RVAL}
        ,
        "normal"={
               test<-Q.Test(X,score)
               STATISTIC<-test$test.statistic
               names(STATISTIC)<-"Q.N"
               PVAL<-test$p.value
               METHOD<-"MARGINAL NORMAL SCORE TEST (assuming the independent component model)"
               RVAL<-list(statistic=STATISTIC,p.value=as.numeric(PVAL),method=METHOD)
               
               RVAL}
        )
        
     PARAMETER<-p
     names(PARAMETER)<-"df" 
    
     if (method=="simulation"){
     Q.simu<-replicate(n.simu,Q.Test(matrix(rnorm(n*p),ncol=p),score)$test.statistic)
     res1$p.value<-mean(Q.simu>res1$statistic)
     PARAMETER<-n.simu
     names(PARAMETER)<-"replications"
     }
     
     if (method=="permutation"){
     Q.simu<-replicate(n.simu,Q.Test(sample(c(1,-1), n, replace = T)*X,score)$test.statistic)
     res1$p.value<-mean(Q.simu>res1$statistic)
     PARAMETER<-n.simu
     names(PARAMETER)<-"replications"
     }
    
    ALTERNATIVE="two.sided"
    NVAL<-paste("c(",paste(mu,collapse=","),")",sep="")
    names(NVAL)<-"location"
    res<-c(res1,list(data.name=DNAME,parameter=PARAMETER,alternative=ALTERNATIVE,null.value=NVAL))
    class(res)<-"htest"
    return(res)
    }

rank.ictest.ics <-
function(X, index = NULL, na.action = na.fail, ...)
    {
    
    DNAME<-deparse(substitute(X))
    if (class(X) != "ics") stop("'ics' must be of class 'ics'") 
    
    Z <- ics.components(X)
    p <- dim(Z)[2]
    if (is.null(index)) index <-  1:p
    DATA<-list(X=Z[,index])
    y <- do.call("rank.ictest", c(DATA, list(...)))
    y$data.name <- DNAME
 
    return(y)    
    }
