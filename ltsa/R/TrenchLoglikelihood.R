`TrenchLoglikelihood` <-
function( r, z )  
{
    EPS<-.Machine$double.eps
    r<-r/r[1]
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
     n<-length(z)
     S<-(out$tr)[1]
     lg<-(out$tr)[2]
    -(n/2)*log(S/n)-(1/2)*lg
}

