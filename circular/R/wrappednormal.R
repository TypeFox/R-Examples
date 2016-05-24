#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   rwrappednormal function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: August, 10, 2006                                  #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.3-3                                           #
#############################################################

rwrappednormal <- function(n, mu=circular(0), rho=NULL, sd=1, control.circular=list()) {
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
   attr(mu, "class") <- attr(mu, "circularp") <- NULL

   if (is.null(rho)) 
      rho <- exp(-sd^2/2)
   if (rho < 0 | rho > 1)
      stop("rho must be between 0 and 1")

   result <- RwrappednormalRad(n, mu, rho)
   result <- conversion.circular(circular(result), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
   return(result)
}

RwrappednormalRad <- function(n, mu, rho) {
   if (rho == 0)
      result <- runif(n, 0, 2*pi)
   else if (rho == 1)
      result <- rep(mu, n)
   else {
      sd <- sqrt(-2 * log(rho))
      result <- rnorm(n, mu, sd) %% (2*pi)
   }
   return(result)
}
   
#############################################################
#                                                           #
#   dwrappednormal function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: May, 31, 2006                                     #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.3-1                                           #
#############################################################

dwrappednormal <- function(x, mu=circular(0), rho=NULL, sd=1, K=NULL, min.k=10) { 
   x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
   mu <- conversion.circular(mu, units="radians", zero=0, rotation="counter")
   attr(x, "class") <- attr(x, "circularp") <- NULL
   attr(mu, "class") <- attr(mu, "circularp") <- NULL

   if (is.null(rho))
      rho <- exp(-sd^2/2)
   if (rho < 0 | rho > 1)
      stop("rho must be between 0 and 1")
   if (length(mu)!=1)
      stop("is implemented only for scalar 'mean'")
   result <- DwrappednormalRad(x, mu, rho, K, min.k)
   return(result)
}

DwrappednormalRad <- function(x, mu, rho, K, min.k=10) {
   var <- -2 * log(rho)
   sd <- sqrt(var)
    
   if (is.null(K)) {
      range <- abs(mu-x)
      K <- (range+6*sqrt(var))%/%(2*pi)+1
      K <- max(min.k, K)
   }
   n <- length(x)
   z <- .Fortran("dwrpnorm",
      as.double(x),
      as.double(mu),
      as.double(sd),
      as.integer(n),
      as.integer(length(mu)),
      as.integer(K),
      d=mat.or.vec(length(mu), n),
      PACKAGE="circular"
   )
   d <- t(z$d/sqrt(var * 2 * pi))
   if (ncol(d)==1)
      d <- c(d)
   return(d)
}

#############################################################
#                                                           #
#   pwrappednormal function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: May, 31, 2006                                     #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2-1                                           #
#############################################################

pwrappednormal <- function(q, mu=circular(0), rho=NULL, sd=1, from=NULL, K=NULL, min.k=10, ...) {
   q <- conversion.circular(q, units="radians", zero=0, rotation="counter", modulo="2pi")
   mu <- conversion.circular(mu, units="radians", zero=0, rotation="counter", modulo="2pi")
   if (is.null(from)) {
      from <- mu - pi
   } else {
      from <- conversion.circular(from, units="radians", zero=0, rotation="counter", modulo="2pi")    
   }

   attr(q, "class") <- attr(q, "circularp") <- NULL
   attr(mu, "class") <- attr(mu, "circularp") <- NULL
   attr(from, "class") <- attr(from, "circularp") <- NULL
   
   n <- length(q)
   if (length(mu) != 1) 
      stop("is implemented only for scalar 'mean'")
   mu <- (mu-from)%%(2*pi)
   q <- (q-from)%%(2*pi)

   if (is.null(rho)) {
      rho <- exp(-sd^2/2)
   }
   if (rho < 0 | rho > 1)
      stop("rho must be between 0 and 1")

   intDwrappednormalRad <- function(q) {
      if (is.na(q)) {    
         return(NA)
      } else {
         return(integrate(DwrappednormalRad, mu=mu, rho=rho, K=K, min.k=min.k, lower=0, upper=q, ...)$value)
      }
   }
   value <- sapply(X=q, FUN=intDwrappednormalRad) 
   return(value)
}

#############################################################
#                                                           #
#   qwrappednormal function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: August, 12, 2010                                  #
#   Copyright (C) 2010 Claudio Agostinelli                  #
#                                                           #
#   Version 0.3-1                                           #
#############################################################

qwrappednormal <- function(p, mu=circular(0), rho=NULL, sd=1, from=NULL, K=NULL, min.k=10, tol=.Machine$double.eps^(0.6), control.circular=list(), ...) {

   epsilon <- 10 * .Machine$double.eps
   if (any(p>1) | any(p<0))
      stop("p must be in [0,1]")
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

   mu <- conversion.circular(mu, units="radians", zero=0, rotation="counter", modulo="2pi")
   if (is.null(from)) {
      from <- mu - pi
   } else {
      from <- conversion.circular(from, units="radians", zero=0, rotation="counter", modulo="2pi")    
   }

   attr(mu, "class") <- attr(mu, "circularp") <- NULL
   attr(from, "class") <- attr(from, "circularp") <- NULL

   n <- length(p)
   if (length(mu) != 1) 
      stop("is implemented only for scalar 'mean'")
   
   mu <- (mu-from)%%(2*pi)
   if (is.null(rho))
      rho <- exp(-sd^2/2)
   if (rho < 0 | rho > 1)
      stop("rho must be between 0 and 1")

   zeroPwrappednormalRad <- function(x, p, mu, rho, K, min.k) {
      if (is.na(x)) {    
         y <- NA
      } else {   
         y <- integrate(DwrappednormalRad, mu=mu, rho=rho, K=K, min.k=min.k, lower=0, upper=x)$value - p
      }
      return(y)
   }

   value <- rep(NA, length(p))
   sem <- options()$show.error.messages
   options(show.error.messages=FALSE)
   for (i in 1:length(p)) {
         res <- try(uniroot(zeroPwrappednormalRad, p=p[i], mu=mu, rho=rho, K=K, min.k=min.k, lower=0, upper=2*pi-epsilon, tol=tol))
         if (is.list(res)) {
             value[i] <- res$root 
         } else if (p[i] < 10*epsilon) {
             value[i] <- 0
         } else if (p[i] > 1-10*epsilon) {
             value[i] <- 2*pi-epsilon
         }
    }
    options(show.error.messages=sem)
    value <- value + from
    value <- conversion.circular(circular(value), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
    return(value)
}
