`mv.1sample.test` <-
function(X, mu=0, score="identity", stand="outer", method = "approximation", n.simu = 1000, na.action=na.fail)
    {
    DNAME<-deparse(substitute(X))
    
    score <- match.arg(score,c("identity","sign","rank"))
    stand <- match.arg(stand,c("inner","outer"))
    method <- match.arg(method,c("approximation","signchange"))
    X<-na.action(X)
    
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric")
    X<-as.matrix(X)
    
    n<-dim(X)[1]
    p<-dim(X)[2]
    if (length(mu)!=p) mu=mu*rep(1,p)
    
    X.centered<-sweep(X,2,mu,"-")
    
      res1<-switch(score,
        "identity"={
               id.loc(X.centered, method=method, n.simu=n.simu)
               }
        ,
        "sign"={
               switch(stand,
                    "outer" = {
                    ssloc.outer(X.centered,method=method,n.simu=n.simu)
                    }
               ,
                    "inner" = {
                    ssloc.inner(X.centered,method=method,n.simu=n.simu)
                    }
                    )
                    }
        ,
        "rank"={
               switch(stand,
                    "outer" = {
                    srloc.outer(X.centered,method=method,n.simu=n.simu)
                    }
               ,
                    "inner" = {
                    srloc.inner(X.centered,method=method,n.simu=n.simu)
                    }
                    )
                    }
        )
        
    NVAL<-paste("c(",paste(mu,collapse=","),")",sep="")
    names(NVAL)<-"location"
    ALTERNATIVE <- "two.sided"
    
    res<-c(res1,list(data.name=DNAME,alternative=ALTERNATIVE,null.value=NVAL))
    class(res) <- "htest"    
    return(res)
    }
