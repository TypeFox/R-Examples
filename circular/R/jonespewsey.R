#############################################################
#                                                           #
#   djonespewsey function                                   #
#   Author: Federico Rotolo                                 #
#   Email: federico.rotolo@stat.unipd.it                    #
#   Date: October, 05, 2010                                 #
#   Copyright (C) 2010 Federico Rotolo                      #
#                                                           #
#   Version                                                 #
#############################################################


djonespewsey <-	function(x, mu=NULL, kappa=NULL, psi=NULL){
  if (is.null(mu) || length(mu)!=1)
    stop("the mean direction parameter 'mu' is mandatory and it must have length 1")
  if (is.null(kappa) || length(kappa)!=1)
    stop("the concentration  parameter 'kappa' is mandatory and it must have length 1")
  if (is.null(psi) || length(psi)!=1)
    stop("the parameter 'psi' is mandatory and it must have length 1")

  if(kappa<0){stop("kappa must be non negative")}

  x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
  mu <- conversion.circular(mu, units="radians", zero=0, rotation="counter")

  mu <- as.vector(mu)
  kappa <- as.vector(kappa)
  psi <- as.vector(psi)
  attr(x, "class") <- attr(x, "circularp") <-  NULL
  attr(mu, "class") <- attr(mu, "circularp") <-  NULL    
  
  DjonespewseyRad(x, mu, kappa, psi)
}


DjonespewseyRad <- function(x, mu, kappa, psi){
	ker<- function(x){ (cosh(kappa*psi)+sinh(kappa*psi)*cos(x-mu))^(1/psi) / (2*pi*cosh(kappa*psi))}
	ncost<-integrate(ker,0,2*pi)$value
	dens<-ker(x)/ncost
	return(dens)
}

