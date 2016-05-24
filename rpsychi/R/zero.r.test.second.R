zero.r.test.second <- function(r, n, sig.level=.05, digits=3){
    r.sem <- (1-r^2)/n
    r.lower <- tanh(atanh(r) +  qnorm(sig.level/2, lower.tail=TRUE)/sqrt(n-3))
    r.upper <- tanh(atanh(r) +  qnorm(sig.level/2, lower.tail=FALSE)/sqrt(n-3))    
    correlation <- c(es = r, lower=r.lower, upper = r.upper, std = r.sem)
    
    
    criterion.power <- c(
            small = power.r(delta=.1, n=n, sig.level=sig.level), 
            medium= power.r(delta=.3, n=n, sig.level=sig.level), 
            large = power.r(delta=.5, n=n, sig.level=sig.level)
    )
        
    output <- list(correlation = correlation, power = criterion.power)
    output <- sapply(output, round, digits)
    return(output)
}
