EW <- function(x,y,beta,tau,h,T0,max.iter,eps) {
  # ---------------------------------------------------------------
  # FUNCTION
  #  Compute expo weight
  # INPUT
  #   x: matrix n,p
  #   y: vector n
  #   beta, tau, h, T0: scalars
  #   max.iter : positive integer number (maximal number of
  #              iterations)
  #   eps  : positive real number (precision)
  # OUTPUT
  #  LEW: p vector
  # CALLED BY
  #   calcModLasso
  # ---------------------------------------------------------------
  n<-length(y)
  p<-dim(x)[2]
  XX<-t(x)%*%x
  Xy<-t(x)%*%y
  a<-2*h*beta/n
  b<-sqrt(2*h)
  c<-2*h
  d<-tau**2
  L<-seq(p)*0
  k<-round(T0/h)
  alea <- rnorm(p*k,mean=0,sd=1)
  #
  # Call Fortran subroutines
  
  veutlw <-  0
  LEW<-seq(p)*0
  
    L <- .Fortran("bouclet",
                  as.integer(p),as.integer(k),
                  as.integer(veutlw),
                  as.double(a),as.double(b),as.double(c),as.double(d),
                  as.double(h),
                  as.double(Xy),as.double(XX),
                  as.double(alea),
                  out=as.double(L),
                  as.double(LEW))$out
  
  #ADD 05/12/2011
  if (all(is.na(L)))
    stop("Error in loops for the method EW: all results are Na")

  crit<-eps
  veutlw <-  1
  i <- 0 #iteration index
  while ( (crit >= eps) && (i  <= max.iter))
    {
    LEWold <- LEW
    alea <- rnorm(p*k,mean=0,sd=1)
      res <- .Fortran("bouclet",
                  as.integer(p),as.integer(k),
                  as.integer(veutlw),
                  as.double(a),as.double(b),as.double(c),as.double(d),
                  as.double(h),
                  as.double(Xy),as.double(XX),
                  as.double(alea),
                  out1=as.double(L),
                  out2=as.double(LEW))
 
    
    L <- res$out1
    LEW <- res$out2
    i <- i + 1
    if (i > 1)	
      crit <- sum((((1-1/i)*LEW-LEWold))**2)/(sum((LEWold)**2)+10**(-10))
  }
    if (i==max.iter) warning("EWmax.iter reached")
  LEW <- LEW/i/T0
  return(LEW)
}
