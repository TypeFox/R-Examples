`exactLoglikelihood` <-
function(r, z, innovationVarianceQ=TRUE)  
{
    EPS<-.Machine$double.eps
    r0 <- r[1] #need to input autocorrelation function to "trenchQFR"
    r <- r/r[1]
    if (r[1] < EPS) 
        stop("error: r[1] the variance is <= 0") 
    if (length(r) != length(z)) 
        stop("error: r and z have unequal lengths") 
     out<-.C("trenchQFR", as.double(r), as.integer(length(r)),
            as.double(z), as.integer(length(z)), EPS,
            tr = array(0, dim=c(1,2 )), fault = as.integer(1),  
      PACKAGE="ltsa" )
     fault<-out$fault
     if (fault == 1) 
        stop("error: matrix is not p.d.")
     n <- length(z)
     S <-(out$tr)[1]
     lg<-(out$tr)[2]
     if (innovationVarianceQ) {
       LL <- -0.5*n*(1+log(2*pi))-0.5*(n*log(S/n)+lg)
       sigmaSq <- S/(n*r0)
     }
     else {
       LL <- sigmaSq <- NA
     }
     list(LL=LL, sigmaSq=sigmaSq)  
}

