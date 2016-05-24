mv.2sample.est <- function(X,g,score="identity", stand="outer", maxiter=100, eps=1e-6, na.action=na.fail, ...)
    {
    dname<-deparse(substitute(X))
    gname<-deparse(substitute(g))
    
    DNAME<-paste(dname,"by",gname)
    
    score <- match.arg(score,c("identity","sign","rank"))
    stand <- match.arg(stand,c("inner","outer"))
    
    if (dim(X)[1]!= length(g)) stop("dimensions of 'X' and 'g' do not match")
    
    COLNAMES.X<-colnames(X)
    
    X<-na.action(data.frame(g=g,X=X))
    
    g<-X$g
    X<-X[,-1]
    
    
    
    if(nlevels(g)!=2) stop("'g' must have two levels")
    
    n.g<-tapply(g,g,length)
    
    if (min(n.g)<2) stop("for each level of 'g' at least two observations must be available")
    
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric")
    
    X<-split.data.frame(X,g)
    
  
    X1<-X[[1]]
    X2<-X[[2]]
    
    
    if (!is.null(COLNAMES.X)){
    colnames(X1)<-COLNAMES.X
    colnames(X2)<-COLNAMES.X
    }
  
    n.X1<-n.g[1]
    n.X2<-n.g[2]

    p<-dim(X1)[2]
    n <- n.X1+n.X2
    
    if(p<2) stop("'X' must be at least bivariate")
    
    #return(list(X1=X1,X2=X2,n=n,p=p))
    
    res1<-switch(score,
        "identity"={
               location<-colMeans(X1)-colMeans(X2)
               covPooled <- ((n.X1-1)/(n.X1+n.X2-2))*cov(X1) + ((n.X2-1)/(n.X1+n.X2-2))*cov(X2)
               scatter<- (1/n.X1+ 1/n.X2)*covPooled
               list(location=location, vcov=scatter, est.name= "difference between sample mean vectors")}
        ,
        "sign"={
               switch(stand,
                    "outer" = {
                    sign.est.outer(X1=X1,X2=X2,p=p,n1=n.X1,n2=n.X2, maxiter=maxiter, eps=eps,...)
                    }
               ,
                    "inner" = {
                    sign.est.inner(X1=as.matrix(X1),X2=as.matrix(X2),p=p,n1=n.X1,n2=n.X2, maxiter=maxiter, eps=eps, ...)
                    }
                    )
                    }
        ,
        "rank"={
               switch(stand,
                    "outer" = {
                    rank.est.outer(X1=as.matrix(X1),X2=as.matrix(X2),p=p,n1=n.X1,n2=n.X2, maxiter=maxiter, eps=eps, ...)
                    }
               ,
                    "inner" = {
                    rank.est.inner(X1=as.matrix(X1),X2=as.matrix(X2),p=p,n1=n.X1,n2=n.X2, maxiter=maxiter, eps=eps,...)
                    }
                    )
                    }
        )
    res1<-c(res1,list(dname=DNAME))
    class(res1) <- "mvloc"    
    return(res1)

    
    
    }
    
