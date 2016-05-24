#############################################################
#                                                           #
#   dgenvonmises function                                   #
#   Author: Federico Rotolo                                 #
#   Email: federico.rotolo@stat.unipd.it                    #
#   Date: October, 05, 2010                                 #
#   Copyright (C) 2010 Federico Rotolo                      #
#                                                           #
#   Version                                                 #
#############################################################

dgenvonmises <- function (x, mu1=NULL, mu2=NULL, kappa1=NULL, kappa2=NULL) {
  if (is.null(mu1) || length(mu1)!=1 || is.null(mu2) || length(mu2)!=1)
    stop("the mean direction parameters 'mu1' and 'mu2' are mandatory and it must have length 1")
  if (is.null(kappa1) || length(kappa1)!=1 || is.null(kappa2) || length(kappa2)!=1)
    stop("the concentration direction parameters 'kappa1' and 'kappa2' are mandatory and it must have length 1")

  if((kappa1<0)||(kappa2<0)){stop("'kappa1' and 'kappa2' must be non negative")}
  if((mu2%/%pi)%%2!=0){stop("'mu2' must be <pi")}

  x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
  mu1 <- conversion.circular(mu1, units="radians", zero=0, rotation="counter")
  mu2 <- conversion.circular(mu2, units="radians", zero=0, rotation="counter")

  mu1 <- as.vector(mu1)
  mu2 <- as.vector(mu2)
  kappa1 <- as.vector(kappa1)
  kappa2 <- as.vector(kappa2)
  attr(x, "class") <- attr(x, "circularp") <-  NULL
  attr(mu1, "class") <- attr(mu1, "circularp") <-  NULL    
  attr(mu2, "class") <- attr(mu2, "circularp") <-  NULL    
  
  DgenvonmisesRad(x, mu1, mu2, kappa1, kappa2)
}


DgenvonmisesRad <-function(x, mu1, mu2, kappa1, kappa2){
	d=(mu1-mu2)%%pi
	num<- exp( kappa1*cos(x-mu1) + kappa2*cos(2*(x-mu2)) )
	den<-integrate(function(x) {
		exp(kappa1*cos(x)+kappa2*cos(2*(x+d)))
	},0,2*pi)$value

	dens<-num/den
	return(dens)
}

