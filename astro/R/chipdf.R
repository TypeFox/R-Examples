chipdf = function(X, k){
    
    p = (1/((2^(k/2))*(gamma(k/2)))) * ((X)^((k/2)-1)) * (exp(-(X)/2))
    
    if(any(X < 0)){p[X<0] = 0}
    
    return(p)
    
}
