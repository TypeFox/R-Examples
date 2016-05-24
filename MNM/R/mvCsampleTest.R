`mv.Csample.test` <-
function(X, g, score="identity", stand="outer", method = "approximation", n.simu = 1000, na.action=na.fail,...)
    {
    DNAME=paste(deparse(substitute(X)),"by",deparse(substitute(g)))
    
    score <- match.arg(score,c("identity","sign","rank"))
    stand <- match.arg(stand,c("inner","outer"))
    method <- match.arg(method,c("approximation","permutation"))
    
    if (length(g)!= dim(X)[1]) stop("'g' must have as many elements as 'X' rows")
    
    DATA <- data.frame(g=g)
    DATA$X <- as.matrix(X)
    DATA<-na.action(DATA)
    
    X <- DATA$X
    g <- DATA$g
    
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric")
    
    p<-dim(X)[2]
    if (p<2) stop("'X' must be at least bivariate")
    n<-dim(X)[1]
    
    if (!is.factor(g)) stop("'g' must be a factor")
    
    n.g<-length(g)
    if (n.g!= n) stop("'g' must have as many elements as 'X' rows")
    
    if (nlevels(g)<2) stop("'g' must have at least two levels")
    
    if (min(by(g,g,length))<2) stop("each level of 'g' must have at least two observations")
    
    
      res1<-switch(score,
        "identity"={
               hot.csample(X,g,method=method,n.simu=n.simu)
               }
        ,
        "sign"={
               switch(stand,
                    "outer" = {
                    CssTestOuter(X,g,method=method,n.simu=n.simu,...)
                    }
               ,
                    "inner" = {
                    CssTestInner(X,g,method=method,n.simu=n.simu,...)
                    }
                    )
                    }
        ,
        "rank"={
               switch(stand,
                    "outer" = {
                    CsrTestOuter(X,g,method=method,n.simu=n.simu,...)
                    }
               ,
                    "inner" = {
                    CsrTestInner(X,g,method=method,n.simu=n.simu,...)
                    }
                    )
                    }
        )
        
    NVAL<-paste("c(",paste(rep(0,p),collapse=","),")",sep="")
    names(NVAL)<-"location difference between some groups"
    ALTERNATIVE <- "two.sided"
    
    res<-c(res1,list(data.name=DNAME,alternative=ALTERNATIVE,null.value=NVAL))
    class(res) <- "htest"    
    return(res)
    }
