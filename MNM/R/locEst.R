`mv.1sample.est` <-
function(X, score="identity", stand="outer", maxiter=100, eps=1e-6, na.action=na.fail, ...)
    {
    dname<-deparse(substitute(X))
    
    score <- match.arg(score,c("identity","sign","rank"))
    stand <- match.arg(stand,c("inner","outer"))
    
    X<-na.action(X)
    
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric")
    X<-as.matrix(X)
    n<-dim(X)[1]
    p<-dim(X)[2]
    
    res1<-switch(score,
        "identity"={
               location<-colMeans(X)
               scatter<- cov(X)/n
               list(location=location, vcov=scatter, est.name= "sample mean vector")}
        ,
        "sign"={
               switch(stand,
                    "outer" = {
                    SIGNS <-spatial.sign(X, center = TRUE, shape = FALSE, maxiter = maxiter, eps = eps, ...)
                    location<-attr(SIGNS,"center")
                    attr(SIGNS,"center")<-NULL
                    attr(SIGNS,"shape")<-NULL
                    r<-SpatialNP:::norm(sweep(X,2,location))
                    w.SIGNS<- SIGNS/sqrt(r)
                    r.sum<-sum(1/r)
                    #A<- (diag(r.sum,p)- t(w.SIGNS) %*% w.SIGNS)/n
                    #B<- t(SIGNS) %*% SIGNS / n
                    A<- (diag(r.sum,p)- crossprod(w.SIGNS))/n
                    B<- crossprod(SIGNS) / n
                    A.inv<-solve(A)
                    scatter <- A.inv %*% B %*% A.inv/n
                    list(location=location, vcov=scatter, est.name= "spatial median")
                    }
               ,
                    "inner" = {
                    SIGNS <-spatial.sign(X, center = TRUE, shape = TRUE, maxiter = maxiter, eps.scale = eps, eps.center=eps, ...)
                    location<-attr(SIGNS,"center")
                    C<-attr(SIGNS,"shape")
                    attr(SIGNS,"center")<-NULL
                    attr(SIGNS,"shape")<-NULL
                    C.inv<-solve(C)
                    H<- mat.sqrt(C.inv)
                    X.inner <- sweep(X,2,location) %*% H
                    r<-SpatialNP:::norm(X.inner)
                    w.SIGNS<- SIGNS/sqrt(r)
                    r.sum<-sum(1/r)
                    #A<- (diag(r.sum,p)- t(w.SIGNS) %*% w.SIGNS)/n
                    #B<- t(SIGNS) %*% SIGNS / n
                    A<- (diag(r.sum,p)- crossprod(w.SIGNS))/n
                    B<- crossprod(SIGNS) / n
                    H.inv <- solve(H)
                    A.inv <- solve(A)
                    scatter <- H.inv %*% A.inv %*% B %*% A.inv %*% H.inv /n
                    list(location=location, vcov=scatter, est.name= "equivariant spatial median")
                    }
                    )
                    }
        ,
        "rank"={
               switch(stand,
                    "outer" = {
                    SIGNRANKS<-spatial.signrank(X, center=TRUE, shape=diag(1,p), maxiter = maxiter, eps = eps, ...)
                    location<-as.vector(attr(SIGNRANKS,"center"))
                    attr(SIGNRANKS,"center")<-NULL
                    attr(SIGNRANKS,"shape")<-NULL
                    Xsums<- pair.sum(sweep(X,2,location))
                    r<-SpatialNP:::norm(Xsums)
                    w.Xsums<- Xsums/(r^1.5)
                    r.sum<-sum(1/r)
                    #A<- (diag(r.sum,p)- t(w.Xsums) %*% w.Xsums)/(length(r))
                    #B<- t(SIGNRANKS) %*% SIGNRANKS / n
                    A<- (diag(r.sum,p)- crossprod(w.Xsums))/(length(r))
                    B<- crossprod(SIGNRANKS) / n
                    A.inv<-solve(A)
                    scatter <- A.inv %*% B %*% A.inv/n
                    res2<-list(location=location, vcov=scatter, est.name= "spatial Hodges-Lehmann estimator")
                    }
               ,
                    "inner" = {
                    SIGNRANKS<-spatial.signrank(X, center=TRUE, shape=TRUE, maxiter = maxiter, eps = eps, ...)
                    location<-attr(SIGNRANKS,"center")
                    C<-attr(SIGNRANKS,"shape")
                    C.inv<-solve(C)
                    attr(SIGNRANKS,"center")<-NULL
                    attr(SIGNRANKS,"shape")<-NULL
                    H<- mat.sqrt(C.inv)
                    X.inner <- sweep(X,2,location) %*% H
                    X.sums<- pair.sum(X.inner)
                    r<-SpatialNP:::norm(X.sums)
                    w.Xsums<- X.sums/(r^1.5)
                    r.sum<-sum(1/r)
                    #A<- (diag(r.sum,p)- t(w.Xsums) %*% w.Xsums)/(length(r))
                    #B<- t(SIGNRANKS) %*% SIGNRANKS / n
                    A<- (diag(r.sum,p)- crossprod(w.Xsums))/(length(r))
                    B<- crossprod(SIGNRANKS) / n
                    A.inv<-solve(A)
                    H.inv <- solve(H)
                    scatter <- H.inv %*% A.inv %*% B %*% A.inv %*% H.inv /n
                    res2 <- list(location=location, vcov=scatter, est.name= "equivariant spatial Hodges-Lehmann estimator")
                    }
                    )
                    }
        )
    res1<-c(res1,list(dname=dname))
    class(res1) <- "mvloc"    
    return(res1)
    }
