chicdf = function(X, k){
    
    #incgamma = pgamma(X/2, k/2, lower=TRUE)*gamma(k/2)
    
    #c  = (1/gamma(k/2)) * incgamma
    
    c = pgamma(X/2, k/2, lower.tail=TRUE)
    
    return(c)
    
}
