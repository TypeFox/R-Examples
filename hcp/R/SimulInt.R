SimulInt <-
function(x.bar, S, a, n, level=0.95){
    ## How to force a condition on the arguments:
    if(length(x.bar) != length(a)){ 
        stop("a and x.bar have to be of equal length")
    }
    p <- length(x.bar) 
    m <- t(a) %*% x.bar
    Sstuff <- (t(a) %*% S %*% a)/n
    Fstuff <- qf(level, p,n-p)*p*(n-1)/(n-p)
    return(c(m - sqrt(Fstuff*Sstuff), m + sqrt(Fstuff*Sstuff)) )
}
