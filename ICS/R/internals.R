### Internal functions for cov4 and ics

### covariance matrix based on 4th moments wrt to the mean vector
### subroutine of cov4
###

.cov4moments.mean<-function(X)
    {
    p<-dim(X)[2]
    n<-dim(X)[1]
    data.centered<-sweep(X,2,colMeans(X),"-")
    Sigma.data.sqrt<-mat.sqrt(cov(X)) 
    radius<-sqrt(rowSums((data.centered %*% solve(Sigma.data.sqrt))^2))
    y<-radius*data.centered
    V<-(1/(n*(p+2)))*crossprod(y) 
    return(V) 
    }

### covariance matrix based on 4th moments wrt to origin
### subroutine of cov4
###

.cov4moments.origin<-function(X)
    {
    p<-dim(X)[2]
    n<-dim(X)[1]
    Sigma.data.sqrt<-mat.sqrt(covOrigin(X)) 
    radius<-sqrt(rowSums((X %*% solve(Sigma.data.sqrt))^2))
    V<-(1/(p+2))*covOrigin(radius*X)  
    return(V) 
    }
    

### Sign of the maximum element of a vector
### returns 1 if the absolute largest value is positive and -1 otherwise
### subroutine in ics
     
.sign.max<-function(x)
 {
 ifelse(identical(max(x),max(abs(x))),1,-1)
 }
        
