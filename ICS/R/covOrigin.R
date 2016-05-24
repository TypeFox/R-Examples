### computes the regular covariance matrix with respect to a fixed location
### the default is the origin
###

`covOrigin` <-
function(X,location=NULL,na.action=na.fail)   
    {
    X<-na.action(X)
    X.matrix<-as.matrix(X)
    
    p<-dim(X.matrix)[2]
    
    if(!is.null(location))
        {
        if(length(location)!=p) stop("'location' is of wrong dimension")
        X.matrix<-sweep(X.matrix,2,location)
        }
    
    n<-dim(X.matrix)[1]
    res<-(1/n)*crossprod(X.matrix)
    return(res)
    }
