pprodnormalMeeker<- function(q, mu.x, mu.y, se.x, se.y, rho=0, lower.tail=TRUE){
    mu.a <- mu.x/se.x
    mu.b <- mu.y/se.y
    q <- q/(se.x*se.y)
    s.a.on.b <- sqrt(1-rho^2) #SD of b conditional on a
    gx=function(x, z) {
        ##Defining the integrand (Meeker & Escobar, 1995, p. 273)
        mu.a.on.b <- mu.a+rho*(x-mu.b) #mean of a conditional on b
        integ <- pnorm(sign(x)*(z / x-mu.a.on.b)/s.a.on.b)*dnorm(x-mu.b)
        return(integ)
    }
    ans<-integrate(gx,lower=-Inf,upper=Inf, z=q )
    percentile <- ans$value
    error <- ans$abs.error
    if (!lower.tail)
      percentile <- 1-percentile
    return(list(p=percentile,error=error))
}
