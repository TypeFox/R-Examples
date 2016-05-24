DAGOSTINO <- function(data){

MOMENTS <- function(data,r) sum((data-mean(data))^r)/length(data)

# Peter
SKEW <- function(data) MOMENTS(data,3)/(MOMENTS(data,2)*sqrt(MOMENTS(data,2)))

# Peter
KURTOSIS <- function(data) MOMENTS(data,4)/(MOMENTS(data,2)*MOMENTS(data,2))

    cat("D'Agostino Test\n")
    n <- length(data)
    cat("    Skewness\n")
    cat("      Skewness coefficient:",sqrtb1 <- SKEW(data),"\n")
    if(n>8){
        y <- sqrtb1*sqrt((n+1)*(n+3)/(6*(n-2)))
        beta2 <- 3*(n*n+27*n-70)*(n+1)*(n+3)/((n-2)*(n+5)*(n+7)*(n+9))
        w <- sqrt(-1+sqrt(2*(beta2-1)))
        delta <- 1/sqrt(log(w))
        ALPHA <- sqrt(2/(w*w-1))
        cat("      Statistics:",zb1 <-
delta*log(y/ALPHA+sqrt((y/ALPHA)^2+1)),"\n")
        cat("      p-value:",2*(1-pnorm(abs(zb1))),"\n")
    }else
        cat("    The skewness test requieres sample size n>8\n")
    cat("    Kurtosis\n")
    cat("      The kurtosis coefficient:",b2 <- KURTOSIS(data),"\n")
    if(n>19){
       meanb2 <- 3*(n-1)/(n+1)
       varb2 <- 24*n*(n-2)*(n-3)/((n+1)*(n+1)*(n+3)*(n+5))
       x <- (b2-meanb2)/sqrt(varb2)
       moment <-
6*(n*n-5*n+2)/((n+7)*(n+9))*sqrt(6*(n+3)*(n+5)/(n*(n-2)*(n-3)))
       a <- 5+8/moment*(2/moment+sqrt(1+4/(moment*moment)))
       cat("      Statistics:",zb2 <-
(1-2/(9*a)-((1-2/a)/(1+x*sqrt(2/(a-4))))^(1/3))/sqrt(2/(9*a)),"\n")
       cat("      p-value:",2*(1-pnorm(abs(zb2))),"\n")
       cat("    Omnibus Test\n")
       cat("      Chi-squared:",k2 <- zb1*zb1+zb2*zb2,"\n")
       cat("      Degree of freedom: 2\n")
       cat("      p-value:",probji2 <- 1-pchisq(k2,2),"\n")
    }else
       cat("    Kurtosis and Omnibus Test require n>=20\n")
}
