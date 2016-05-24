Hypergeometric2F1 <- function(a,b,c,x)
  {

    ## Method 3
    ## numerical integration of results on pages 231 and 233 of Mathai 1999
    ## only used when b = 1/2, c=3/2
    ## or when b = 1, c = 5/2
    ret <- NA
    if(is.nan(x)) {
      ret <- 0
    } else if(x==0) ret <- 1
    if(is.na(ret))
      {
        ## this seems to be the best, in terms of speed and accuracy
        if(b==1 && c==5/2) {
          ## use Mathai 1999, page 231, lemma 2.4.1
          integrand1 <- function(u) u^(-2*a+1)*sqrt(u^2-(1-x))
          result <- integrate(integrand1, lower=sqrt(1-x),upper=1,rel.tol=1e-12)$value
          result <- result*3/(x*sqrt(x))
          ret <- result
        } else if(b==1/2 && c==3/2)
          {
            if(x==1) {
              ret <- sqrt(pi)*gamma(1-a)/(2*gamma(3/2-a))
            } else {
              ## use Mathai 1999, page 233, lemma 2.4.2
              integrand2 <- function(u) u^(2*a-2)*asin(sqrt(1-x)/u)
              result <- integrate(integrand2, lower=sqrt(1-x),upper=1,rel.tol=1e-12)$value
              
              result <- result*(2*a-1)/sqrt(1-x)^(2*a-1)+pi/2-1/sqrt(1-x)^(2*a-1)*asin(sqrt(1-x))
              result <- result/sqrt(x)
              ret <- result
            }
          }
      }
    ##cat("hypgeo =",ret,"\n")
    ##if(is.nan(ret)) cat(a,b,c,x,"error\n")
    if(is.nan(ret)) ret <- 0  ## if this happens 2F1 will be multiplied by 0 anyway
    ret
}
