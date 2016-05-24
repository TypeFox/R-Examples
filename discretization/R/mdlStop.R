mdlStop <-
function(ci,y,entropy){
    n <- length(y)
    es <- ent(y)
    left <- 1:ci; right <- (ci+1):n     
    gain <- es-entropy   
    l0 <- levels(factor(y))
    l1 <- levels(factor(y[left])); l2 <- levels(factor(y[right]))         
    k <- length(l0)
    k1 <- length(l1); k2 <- length(l2)
    delta <- mylog(3^k-2)-(k*es-k1*ent(y[left])-k2*ent(y[right]))
    cond <- mylog(n-1)/n+delta/n        
    if(gain<cond) return (NULL)
    return(gain)
}
