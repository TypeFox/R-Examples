### Function for the two sample Hodges Lehmann estimator
### -> internal use only
###

HL.est <- function(X1,X2,...)
    {
    n1<-dim(X1)[1]
    n2<-dim(X2)[1]
    
    index1<- rep(1:n1,n2)
    index2<- rep(1:n2,each=n1)
    
    HL <- spatial.median(X1[index1,]-X2[index2,],...)
    return(HL)
    }
    
# Alternative function
    
HL.est2 <- function(X1,X2,...)
    {
    m <- nrow(X1)*nrow(X2)
    Z <- matrix( rep(t(X1), nrow(X2)), byrow=TRUE, nrow=m) - matrix( rep(X2, each=nrow(X1),byrow=TRUE),nrow=m) 
    HL <- spatial.median(Z,...)
    return(HL)
    }

### Function for the joint estimation of the transformation matrix X
### in the two sample case when spatial signs are used
###

Hsign<-function(X1, X2, eps = 1e-06, maxiter = 500, ...)
    {
    p<-dim(X1)[2]
    n1<-dim(X1)[1]
    n2<-dim(X2)[1]
    C<-diag(p)
    
    differ<-Inf
    iter<-0
    while (differ>eps)
        {
        if (maxiter<= iter) stop("'maxiter' reached without convergence") 
        C.sqrt <- mat.sqrt(C)
        C.inv <- solve(C)
        C.inv.sqrt <- mat.sqrt(C.inv)
    
        Y1 <- X1 %*% C.inv.sqrt
        Y2 <- X2 %*% C.inv.sqrt
    
        smed.X1<-spatial.median(Y1, maxiter=maxiter, eps=eps, ...) %*% C.sqrt
        smed.X2<-spatial.median(Y2, maxiter=maxiter, eps=eps, ...) %*% C.sqrt
    
        residuals.X1 <- sweep(X1,2,as.vector(smed.X1)) %*% C.inv.sqrt
        residuals.X2 <- sweep(X2,2,as.vector(smed.X2)) %*% C.inv.sqrt
    
        ss.r.X1 <- spatial.sign(residuals.X1, FALSE, FALSE)
        ss.r.X2 <- spatial.sign(residuals.X2, FALSE, FALSE)
    
        B1<- t(ss.r.X1) %*% ss.r.X1
        B2<- t(ss.r.X2) %*% ss.r.X2
    
        B<-(1/(n1+n2))*(B1+B2)
    
        C.new <- p * C.sqrt %*% B %*% C.sqrt
        differ<-frobenius.norm(C.new-C)
        C<-C.new
        iter<-iter+1
        }
    return(C)
    }

### Function for the joint estimation of the transformation matrix X
### in the two sample case when spatial ranks are used
###

Hrank<-function(X1, X2, eps = 1e-06, maxiter = 500, ...)
    {
    p<-dim(X1)[2]
    n1<-dim(X1)[1]
    n2<-dim(X2)[1]
    C.mat<-diag(p)
    
    differ<-Inf
    iter<-0
    while (differ>eps)
        {
        if (maxiter<= iter) stop("'maxiter' reached without convergence") 
        C.sqrt <- mat.sqrt(C.mat)
        C.inv.sqrt<-solve(C.sqrt)
    
        Y1 <- X1 %*% C.inv.sqrt
        Y2 <- X2 %*% C.inv.sqrt
        
        delta<-HL.est(Y1,Y2,maxiter=maxiter, eps=eps, ...) %*% C.sqrt
        
        Y2b <- sweep(X2,2,delta,"+") %*% C.inv.sqrt
    
        Ycomb <- rbind(Y1,Y2b)
        
        RANKS <- spatial.rank(Ycomb, shape = FALSE)
       
        
        B <- crossprod(RANKS) /(n1+n2)
    
        
    
        C.new <- C.sqrt %*% B %*% C.sqrt
        C.new <- p * 1/sum(diag(C.new)) * C.new
        
        differ<-frobenius.norm(C.new-C.mat)
        C.mat<-C.new
        iter<-iter+1
        }
    return(C.mat)
    
    }

sign.est.outer<-function(X1,X2, p, n1, n2, maxiter, eps, ...){
                    n<-n1+n2
                    SIGNS.X1 <-spatial.sign(X1, center = TRUE, shape = FALSE, maxiter=maxiter, eps=eps,...)
                    SIGNS.X2 <-spatial.sign(X2, center = TRUE, shape = FALSE, maxiter=maxiter, eps=eps,...)
                    center.X1<-attr(SIGNS.X1,"center")
                    center.X2<-attr(SIGNS.X2,"center")
                    location<- center.X1-center.X2
                    attr(SIGNS.X1,"center")<-NULL
                    attr(SIGNS.X1,"shape")<-NULL
                    attr(SIGNS.X2,"center")<-NULL
                    attr(SIGNS.X2,"shape")<-NULL
                    SIGNS<-rbind(SIGNS.X1,SIGNS.X2)
                    X1.centered <- sweep(X1,2,center.X1)
                    X2.centered <- sweep(X2,2,center.X2)
                    r<-SpatialNP:::norm(as.matrix(rbind(X1.centered,X2.centered)))
                    w.SIGNS<- SIGNS/sqrt(r)
                    r.sum<-sum(1/r)
                    A<- (diag(r.sum,p)- t(w.SIGNS) %*% w.SIGNS)/n
                    B.X1<- t(SIGNS.X1) %*% SIGNS.X1 
                    B.X2<- t(SIGNS.X2) %*% SIGNS.X2
                    B<- (B.X1 + B.X2)/n
                    A.inv<-solve(A)
                    scatter <- (1/n1 + 1/n2)* A.inv %*% B %*% A.inv
                    list(location=location, vcov=scatter, est.name= "difference between spatial medians")
                    }

sign.est.inner<-function(X1, X2, p, n1, n2, maxiter, eps,...){
                    n<-n1+n2
                    C<-Hsign(X1, X2, maxiter=maxiter, eps=eps)
                    SIGNS.X1 <-spatial.sign(X1, center = TRUE, shape = C, maxiter=maxiter, eps=eps,...)
                    SIGNS.X2 <-spatial.sign(X2, center = TRUE, shape = C, maxiter=maxiter, eps=eps,...)
                    center.X1<-attr(SIGNS.X1,"center")
                    center.X2<-attr(SIGNS.X2,"center")
                    location<- center.X1-center.X2
                    C.inv<-solve(C)
                    H<- mat.sqrt(C.inv)
                    X1.inner <- sweep(X1,2,center.X1) %*% H
                    X2.inner <- sweep(X2,2,center.X2) %*% H
                    attr(SIGNS.X1,"center")<-NULL
                    attr(SIGNS.X1,"shape")<-NULL
                    attr(SIGNS.X2,"center")<-NULL
                    attr(SIGNS.X2,"shape")<-NULL
                    SIGNS<-rbind(SIGNS.X1,SIGNS.X2)
                    r<-SpatialNP:::norm(as.matrix(rbind(X1.inner,X2.inner)))
                    w.SIGNS<- SIGNS/sqrt(r)
                    r.sum<-sum(1/r)
                    A<- (diag(r.sum,p)- t(w.SIGNS) %*% w.SIGNS)/n
                    B.X1<- t(SIGNS.X1) %*% SIGNS.X1 
                    B.X2<- t(SIGNS.X2) %*% SIGNS.X2
                    B<- (B.X1 + B.X2)/n
                    H.inv <- solve(H)
                    A.inv <- solve(A)
                    scatter <- (1/n1 + 1/n2)*H.inv %*% A.inv %*% B %*% A.inv %*% H.inv 
                    list(location=location, vcov=scatter, est.name= "difference between equivariant spatial medians")
                    }


rank.est.outer<-function(X1,X2,p,n1,n2,maxiter, eps, ...)              
                    {
                    n <- n1+n2
                    #RANKS.X1 <- spatial.rank(X1, shape=FALSE)
                    #RANKS.X2 <- spatial.rank(X2, shape=FALSE)
                    index1 <- rep(1:n1,n2)
                    index2 <- rep(1:n2,each=n1)     
                    location <- spatial.median(X1[index1,]-X2[index2,], maxiter=maxiter, eps=eps, ...)
                    X2.equalX1 <- sweep(X2,2,location,"+")
                    X1minusX2 <- X1[index1,]-X2.equalX1[index2,]
                    r<-SpatialNP:::norm(X1minusX2)
                    w.X1minusX2<- X1minusX2/(r^1.5)
                    r.sum<-sum(1/r)
                    A<- (diag(r.sum,p)- t(w.X1minusX2) %*% w.X1minusX2)/(length(r))
                    #B1<- t(RANKS.X1) %*% RANKS.X1 
                    #B2<- t(RANKS.X2) %*% RANKS.X2 
                    X1X2<- rbind(X1,X2.equalX1)
                    RANKS.X1X2 <- spatial.rank(X1X2, shape=FALSE)
                    B<- (crossprod(RANKS.X1X2))/n
                    A.inv<-solve(A)
                    scatter <- (1/n1 + 1/n2)*(A.inv %*% B %*% A.inv)
                    list(location=location, vcov=scatter, est.name= "spatial Hodges-Lehmann estimator for location difference")
                    }

rank.est.inner<-function(X1, X2, p, n1, n2, maxiter, eps, ...)
                    {
                    n<-n1+n2
                    C.mat <- Hrank(X1, X2, maxiter=maxiter, eps=eps)
                    C.inv<-solve(C.mat)
                    H<- mat.sqrt(C.inv)
                    H.inv <- solve(H)
                    
                    
                    #RANKS.X1<-spatial.rank(X1, shape=C.amt)
                    #RANKS.X2<-spatial.rank(X2, shape=C.mat)
                    Y1 <- X1 %*% H
                    Y2 <- X2 %*% H
                    index1 <- rep(1:n1,n2)
                    index2 <- rep(1:n2,each=n1)     
                    
                    locationY<- HL.est(Y1,Y2,maxiter=maxiter, eps=eps,...)
                    location <- locationY %*% H.inv
                    
                    Y2.equalY1 <- sweep(Y2,2,locationY,"+")
                    Y1minusY2 <- Y1[index1,]-Y2.equalY1[index2,]
                  
                   
                    r<-SpatialNP:::norm(Y1minusY2)
                    w.Ysums<- Y1minusY2/(r^1.5)
                    r.sum<-sum(1/r)
                    A<- (diag(r.sum,p)- t(w.Ysums) %*% w.Ysums)/(length(r))
                    #B1<- t(RANKS.X1) %*% RANKS.X1 
                    #B2<- t(RANKS.X2) %*% RANKS.X2
                    
                    Y1Y2<- rbind(Y1,Y2.equalY1)
                    RANKS.Y1Y2 <- spatial.rank(Y1Y2, shape=FALSE)
                    B<- crossprod(RANKS.Y1Y2)/n
                    A.inv<-solve(A)
                    scatter <- (1/n1 + 1/n2)*(H.inv %*% A.inv %*% B %*% A.inv %*% H.inv)
                    list(location=location, vcov=scatter, est.name= "equivariant spatial Hodges-Lehmann estimator for location difference")
                    }
