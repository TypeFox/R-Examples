#############################################################
#                                                           #
#   rkatojones function                                     #
#   Author: Federico Rotolo,                                #
#           original code from Kato, S. and Jones, M.C.     #
#   Email: federico.rotolo@stat.unipd.it                    #
#   Date: October, 23, 2010                                 #
#   Copyright (C) 2010 Federico Rotolo                      #
#                                                           #
#   Version                                                 #
#############################################################

rkatojones <- function(n, mu=NULL, nu=NULL, r=NULL, kappa=NULL, control.circular=list()) {
  if (is.null(mu) || length(mu)!=1)
    stop("the mean direction parameter 'mu' is mandatory and it must have length 1")
  if (is.null(nu) || length(nu)!=1)
    stop("the parameter 'nu' is mandatory and it must have length 1")
  if (is.null(r) || length(r)!=1)
    stop("the parameter 'r' is mandatory and it must have length 1")
  if (is.null(kappa) || length(kappa)!=1)
    stop("the parameter 'kappa' is mandatory and it must have length 1")
  if((r<0)||(r>=1)){stop("'r' must be in [0,1)")}
  if(kappa<0){stop("'kappa' must be not negative")}

  if (is.circular(mu)) {
    datacircularp <- circularp(mu)
  } else {
    datacircularp <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")
  }
  dc <- control.circular
  if (is.null(dc$type))
    dc$type <- datacircularp$type
  if (is.null(dc$units))
    dc$units <- datacircularp$units
  if (is.null(dc$template))
    dc$template <- datacircularp$template
  if (is.null(dc$modulo))
    dc$modulo <- datacircularp$modulo
  if (is.null(dc$zero))
    dc$zero <- datacircularp$zero
  if (is.null(dc$rotation))
    dc$rotation <- datacircularp$rotation
   
  mu <- conversion.circular(mu, units="radians", zero=0, rotation="counter")
  mu <- as.vector(mu)
  nu <- conversion.circular(nu, units="radians", zero=0, rotation="counter")
  nu <- as.vector(nu)
  r <- as.vector(r)
  kappa <- as.vector(kappa)
  attr(mu, "class") <- attr(mu, "circularp") <-  NULL
  attr(nu, "class") <- attr(nu, "circularp") <-  NULL
  vm <- rkatojonesRad(n, mu, nu, r, kappa)
  vm <- conversion.circular(circular(vm), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
  return(vm)
}

rkatojonesRad <- function(n, mu, nu, r, kappa) {
   x <- vector(length = n)
   if (kappa) {
	t <- NULL
        j <- 1
        a <- 1 + sqrt(1 + 4 * kappa^2)
        b <- (a - sqrt(2 * a)) / (2 * kappa)
        zeta <- (1 + b^2) / (2 * b)
        while (j <= n) {
           u <- runif(2)
           z <- cos(pi * u[1])
           f <- (zeta * z + 1) / (zeta + z)
           con <- kappa * (zeta - f)
           if ((con * (2 - con) - u[2] > 0) || (log(con / u[2]) + 1 - con >= 0)) {
              u3 <- runif(1)
              t[j] <- sign(u3 - 0.5) * acos(f)
              j <- j+1
           }
        }
   } else {
     t <- runif(n, min = -pi, max = pi) 
   }
   vm <- mu + nu + 2 * atan((1-r) / (1+r) * tan((t-nu) / 2))
   vm <- Arg(exp((1i) * vm))

   return(vm)
}

#############################################################
#                                                           #
#   dkatojones function                                     #
#   Author: Federico Rotolo                                 #
#   Email: federico.rotolo@stat.unipd.it                    #
#   Date: October, 05, 2010                                 #
#   Copyright (C) 2010 Federico Rotolo                      #
#                                                           #
#   Version                                                 #
#############################################################


dkatojones <-function(x, mu=NULL, nu=NULL, r=NULL, kappa=NULL){
  if (is.null(mu) || length(mu)!=1)
    stop("the mean direction parameter 'mu' is mandatory and it must have length 1")
  if (is.null(nu) || length(nu)!=1)
    stop("the parameter 'nu' is mandatory and it must have length 1")
  if (is.null(r) || length(r)!=1)
    stop("the parameter 'r' is mandatory and it must have length 1")
  if (is.null(kappa) || length(kappa)!=1)
    stop("the parameter 'kappa' is mandatory and it must have length 1")
  if((r<0)||(r>=1)){stop("'r' must be in [0,1)")}
  if(kappa<0){stop("'kappa' must be not negative")}

  x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
  mu <- conversion.circular(mu, units="radians", zero=0, rotation="counter")
  nu <- conversion.circular(nu, units="radians", zero=0, rotation="counter")

  mu <- as.vector(mu)
  nu <- as.vector(nu)
  kappa <- as.vector(kappa)
  r <- as.vector(r)
  attr(x, "class") <- attr(x, "circularp") <-  NULL
  attr(mu, "class") <- attr(mu, "circularp") <-  NULL    
  attr(nu, "class") <- attr(nu, "circularp") <-  NULL    
  
  DkatojonesRad(x, mu, nu, r, kappa)
}


DkatojonesRad <-function(x, mu, nu, r, kappa){
	gamma<-mu+nu
	den <- 2*pi*besselI(kappa,0) * (1+r^2-2*r*cos(x-gamma))
	xi<-(r^4+2*r^2*cos(2*nu)+1)^.5
	eta<-mu+Arg(r^2*cos(2*nu)+r^2*sin(2*nu)*1i+1)
	num<-(1-r^2)*exp((kappa*(xi*cos(x-eta)-2*r*cos(nu)))/(1+r^2-2*r*cos(x-gamma)))
	return(num/den)
}
