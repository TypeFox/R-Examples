effectiveAlpha <- function(rho, alpha, N, df)
{
    sigma <- matrix(c(1,rho,rho,1), ncol =2 , nrow = 2)
    q <- abs(qt(1-alpha/2,df))
    numerator <- as.numeric(pmvnorm(lower = c(-Inf, -q), upper=c(-q, q), mean = c(0,0), sigma = sigma, 
                           algorithm=GenzBretz(abseps=.Machine$double.eps)))
    numerator <- numerator +  as.numeric(pmvnorm(lower = c(q, -q), upper=c(Inf, q), mean = c(0,0), sigma = sigma, 
                                        algorithm=GenzBretz(abseps=.Machine$double.eps)))
    denominator <- pnorm(q) - pnorm(-q) 
    out <- numerator/denominator
    return(1-(1-out)^(N-1) * denominator)
}