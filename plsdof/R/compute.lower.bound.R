compute.lower.bound=function(X){
    S=cor(X)
    lower.bound=-1
    lambda<-eigen(S)$values
    bound=FALSE
    if (2*max(lambda)<=sum(lambda)){
        bound=TRUE
        lower.bound<-1+ sum(lambda)/max(lambda)
    } 
    return(list(bound=bound,lower.bound=lower.bound))
}
