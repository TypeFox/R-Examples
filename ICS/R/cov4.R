### Wrapper for the covariance matrix of fourth moments
### calls depending on the location argument different subroutines
###

`cov4` <-
function(X, location="Mean", na.action=na.fail)
    {
    X<-na.action(X)
    X.matrix<-as.matrix(X)
    
    
    p <- dim(X)[2]                                                
    if (p<2) stop("'X' must be at least bivariate")  
    
    if(is.numeric(location))
        {
        if(length(location)!=p) stop("'location' is of wrong dimension")
        X.matrix<-sweep(X.matrix,2,location)
        location="Origin"
        }
    
    loc<-match.arg(location,c("Mean","Origin"))
    if (loc=="Mean")
        V<-.cov4moments.mean(X.matrix)
    if (loc=="Origin")
        V<-.cov4moments.origin(X.matrix)
    return(V)
    }
