acvfPLS<-function(alpha,maxlag){
    rasy <- (2^(-1 + alpha)*pi^(1/2 - alpha)*alpha*gamma(alpha/2))/gamma(1/2 - alpha/2)
    c(1,rasy*(1:maxlag)^(-alpha))
    }
