### covAxis computes a scatter functional that is needed for PAA
### The matrix can be seen as a one step Tyler functional
###

`covAxis` <-
function(X, na.action = na.fail)
    {
    X<-na.action(X)
    X<-as.matrix(X)
    
    n <- dim(X)[1]
    p <- dim(X)[2]                                                
    if (p<2) stop("'X' must be at least bivariate")  
    
    Xmeans <- colMeans(X)
    di<-sqrt(mahalanobis(X,Xmeans,cov(X)))
    X.centered <- sweep(X, 2, Xmeans)
    Y<-sweep(X.centered,1,di,FUN="/")

    v.tilde <- p*crossprod(Y) / n
    return(v.tilde)
    }
