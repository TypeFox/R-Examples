.Q.simu.sign<-function(X,p1,p2,n)
    {
    Z1.signs <- apply(as.matrix(X[,1:p1]),2,sign)
    Z1.ranks <- apply(abs(as.matrix(X[,1:p1])),2,rank)
    
    Z2.signs <- apply(as.matrix(X[,-(1:p1)]),2,sign)
    Z2.ranks <- apply(abs(as.matrix(X[,-(1:p1)])),2,rank)
    
    C.stat <- t(Z1.signs) %*% Z2.signs/n
    n*frobenius.norm(as.matrix(C.stat))^2
    }

.Q.simu.rank<-function(X,p1,p2,n)
    {
    Z1.signs <- apply(as.matrix(X[,1:p1]),2,sign)
    Z1.ranks <- apply(abs(as.matrix(X[,1:p1])),2,rank)
    
    Z2.signs <- apply(as.matrix(X[,-(1:p1)]),2,sign)
    Z2.ranks <- apply(abs(as.matrix(X[,-(1:p1)])),2,rank)
    
    C.stat <- t(Z1.signs*Z1.ranks)%*%(Z2.signs*Z2.ranks) * 2/(n*(n+1)^2)
    n*frobenius.norm(as.matrix(C.stat))^2
    }

.Q.simu.normal<-function(X,p1,p2,n)
    {
    Z1.signs <- apply(as.matrix(X[,1:p1]),2,sign)
    Z1.ranks <- apply(abs(as.matrix(X[,1:p1])),2,rank)
    
    Z2.signs <- apply(as.matrix(X[,-(1:p1)]),2,sign)
    Z2.ranks <- apply(abs(as.matrix(X[,-(1:p1)])),2,rank)
    
    C.stat <- t(Z1.signs*apply((Z1.ranks/(n+1)+1)/2,2,qnorm))%*%(Z2.signs*apply((Z2.ranks/(n+1)+1)/2,2,qnorm))/n
    n*frobenius.norm(as.matrix(C.stat))^2
    }


.Q.perm.sign<-function(Z1.signs,Z1.ranks,Z2.signs,Z2.ranks,n)
    {
    perm.index <- sample(1:n,n)
    Z2.signs <- Z2.signs[perm.index,]
    Z2.ranks <- Z2.ranks[perm.index,]
    
    C.stat <- t(Z1.signs) %*% Z2.signs/n
    n*frobenius.norm(as.matrix(C.stat))^2
    }

.Q.perm.rank<-function(Z1.signs,Z1.ranks,Z2.signs,Z2.ranks,n)
    {
    perm.index <- sample(1:n,n)
    Z2.signs <- Z2.signs[perm.index,]
    Z2.ranks <- Z2.ranks[perm.index,]
    
    C.stat <- t(Z1.signs*Z1.ranks)%*%(Z2.signs*Z2.ranks) * 2/(n*(n+1)^2)
    n*frobenius.norm(as.matrix(C.stat))^2
    }

.Q.perm.normal<-function(Z1.signs,Z1.ranks,Z2.signs,Z2.ranks,n)
    {
    perm.index <- sample(1:n,n)
    Z2.signs <- Z2.signs[perm.index,]
    Z2.ranks <- Z2.ranks[perm.index,]
    
    C.stat <- t(Z1.signs*apply((Z1.ranks/(n+1)+1)/2,2,qnorm))%*%(Z2.signs*apply((Z2.ranks/(n+1)+1)/2,2,qnorm))/n
    n*frobenius.norm(as.matrix(C.stat))^2
    }
