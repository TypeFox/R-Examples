#############################################################
#                                                           #
#   dcarthwrite function                                    #
#   Author: Federico Rotolo                                 #
#   Email: federico.rotolo@stat.unipd.it                    #
#   Date: October, 05, 2010                                 #
#   Copyright (C) 2010 Federico Rotolo                      #
#                                                           #
#   Version                                                 #
#############################################################

dcarthwrite <- function (x, mu=NULL, psi=NULL) {
  if (is.null(mu) || length(mu)!=1)
    stop("the mean direction parameter 'mu' is mandatory and it must have length 1")

  if (is.null(psi) || length(psi)!=1)
    stop("the parameter 'psi' is mandatory and it must have length 1")

  if(psi<0)
    stop("the parameter 'psi' must be non negative")

  x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
  mu <- conversion.circular(mu, units="radians", zero=0, rotation="counter")
  mu <- as.vector(mu)
  psi <- as.vector(psi)

  attr(x, "class") <- attr(x, "circularp") <-  NULL
  attr(mu, "class") <- attr(mu, "circularp") <-  NULL    
  
  DcarthwrightRad(x, mu, psi)
}


DcarthwrightRad <- function(x, mu, psi) {
  cpc<-2^(1/psi-1) * (gamma(1+1/psi))^2 * (1+cos(x-mu))^(1/psi)
  cpc<-cpc/(pi*gamma(1+2/psi))  
  return(cpc)
}

