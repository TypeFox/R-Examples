#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   rvonmises function                                      #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: November, 06, 2013                                #
#   Copyright (C) 2013 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2-6                                           #
#############################################################

rvonmises <- function(n, mu, kappa, control.circular=list()) {
  if (missing(mu) || length(mu)!=1)
    stop("the mean direction parameter 'mu' is mandatory and it must have length 1")
  if (missing(kappa) || length(kappa)!=1)
    stop("the concentration parameter 'kappa' is mandatory and it must have length 1")   
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
  kappa <- as.vector(kappa)
  if (kappa < 0)
    stop("the concentration parameter 'kappa' must be non negative")
  attr(mu, "class") <- attr(mu, "circularp") <-  NULL
  vm <- RvonmisesRad(n, mu, kappa)
  vm <- conversion.circular(circular(vm), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
  return(vm)
}

#RvonmisesRad <- function(n, mu, kappa) {
#   vm <- 1:n
#   a <- 1 + (1 + 4 * (kappa^2))^0.5
#   b <- (a - (2 * a)^0.5)/(2 * kappa)
#   r <- (1 + b^2)/(2 * b)
#  obs <- 1
#   while (obs <= n) {
#      U1 <- runif(1, 0, 1)
#      z <- cos(pi * U1)
#      f <- (1 + r * z)/(r + z)
#      c <- kappa * (r - f)
#      U2 <- runif(1, 0, 1)
#      if (c * (2 - c) - U2 > 0) {
#         U3 <- runif(1, 0, 1)
#         vm[obs] <- sign(U3 - 0.5) * acos(f) + mu
#         vm[obs] <- vm[obs] %% (2 * pi)
#         obs <- obs + 1
#      } else {
#         if (log(c/U2) + 1 - c >= 0) {
#            U3 <- runif(1, 0, 1)
#            vm[obs] <- sign(U3 - 0.5) * acos(f) + mu
#           vm[obs] <- vm[obs] %% (2 * pi)
#           obs <- obs + 1
#         }
#      }
#   }
#   return(vm)
#}

RvonmisesRad <- function(n, mu, kappa) {
   x <- vector(length = n)
   if (kappa) {
     vm <- .C("rvm",
       as.double(x),
       as.integer(n),
       as.double(mu),
       as.double(kappa),
       PACKAGE="circular")[[1]] %% (2 * pi)
   } else {
     vm <- stats::runif(n=n, min=0, max=2*pi)
   }
   return(vm)
}

#############################################################
#                                                           #
#   dvonmises function                                      #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: February, 07, 2013                                #
#   Copyright (C) 2013 Claudio Agostinelli                  #
#                                                           #
#   Version 0.3                                             #
#############################################################

dvonmises <- function (x, mu, kappa, log=FALSE) {
  if (missing(mu) || length(mu)!=1)
    stop("the mean direction parameter 'mu' is mandatory and it must have length 1")
  if (missing(kappa) || length(kappa)!=1)
    stop("the concentration parameter 'kappa' is mandatory and it must have length 1")   
  if (!is.logical(log))
    stop("'log' must be logical")
  x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
  mu <- conversion.circular(mu, units="radians", zero=0, rotation="counter")
  mu <- as.vector(mu)
  kappa <- as.vector(kappa)
  if (kappa < 0)
    stop("the concentration parameter 'kappa' must be non negative")
  attr(x, "class") <- attr(x, "circularp") <-  NULL
  attr(mu, "class") <- attr(mu, "circularp") <-  NULL    
  
  DvonmisesRad(x, mu, kappa, log)
}

DvonmisesRad <- function(x, mu, kappa, log=FALSE) {
  if (log) {
    if (kappa == 0)
      vm <- log(rep(1/(2*pi), length(x)))
    else if (kappa < 100000)
      vm <- -(log(2*pi)+log(besselI(kappa, nu = 0, expon.scaled=TRUE)) + kappa) + kappa*(cos(x - mu))
    else
      vm <- ifelse(((x-mu)%%(2*pi))==0, Inf, -Inf)
  } else {
    if (kappa == 0)
      vm <- rep(1/(2*pi), length(x))
    else if (kappa < 100000)
      vm <- 1/(2 * pi * besselI(x = kappa, nu = 0, expon.scaled = TRUE)) * (exp(cos(x - mu) -1))^kappa
    else
      vm <- ifelse(((x-mu)%%(2*pi))==0, Inf, 0)
  }
  return(vm)
}

#############################################################
#                                                           #
#   pvonmises function                                      #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: August, 12, 2010                                  #
#   Copyright (C) 2010 Claudio Agostinelli                  #
#                                                           #
#   Version 0.4                                             #
#############################################################

pvonmises <- function(q, mu, kappa, from=NULL, tol = 1e-020) {
  if (missing(mu) || length(mu)!=1)
    stop("the mean direction parameter 'mu' is mandatory and it must have length 1")
  if (missing(kappa) || length(kappa)!=1)
    stop("the concentration parameter 'kappa' is mandatory and it must have length 1")   

  q <- conversion.circular(q, units="radians", zero=0, rotation="counter")
  mu <- conversion.circular(mu, units="radians", zero=0, rotation="counter")
  mu <- as.vector(mu)
  kappa <- as.vector(kappa)
  if (kappa < 0)
    stop("the concentration parameter 'kappa' must be non negative")
   if (is.null(from)) {
      from <- mu - pi
   } else {
      from <- conversion.circular(from, units="radians", zero=0, rotation="counter", modulo="2pi")    
   }  
  attr(q, "class") <- attr(q, "circularp") <-  NULL    
  attr(mu, "class") <- attr(mu, "circularp") <-  NULL
  attr(from, "class") <- attr(from, "circularp") <- NULL
  mu <- (mu-from)%%(2*pi)
  q <- (q-from)%%(2*pi)
  
  PvonmisesRad(q, mu, kappa, tol)
}

PvonmisesRad <- function(q, mu, kappa, tol) {    
   q <- q %% (2 * pi)
   n <- length(q)
   mu <- mu %% (2 * pi)
   pvm.mu0 <- function(q, kappa, tol) {
      flag <- TRUE
      p <- 1
      sum <- 0
      while (flag) {
         term <- (besselI(x=kappa, nu=p, expon.scaled = FALSE) * sin(p * q))/p
         sum <- sum + term
         p <- p + 1
         if (abs(term) < tol)
            flag <- FALSE
      }
      return(q/(2 * pi) + sum/(pi * besselI(x=kappa, nu=0, expon.scaled = FALSE)))
   }

   result <- rep(NA, n)
   if (mu == 0) {
      for (i in 1:n) {
         result[i] <- pvm.mu0(q[i], kappa, tol)
      }
   } else {
      for (i in 1:n) {   
         if (q[i] <= mu) {
            upper <- (q[i] - mu) %% (2 * pi)
            if (upper == 0)
               upper <- 2 * pi
            lower <- ( - mu) %% (2 * pi)
            result[i] <- pvm.mu0(upper, kappa, tol) - pvm.mu0(lower, kappa, tol)
         } else {
            upper <- q[i] - mu
            lower <- mu %% (2 * pi)
            result[i] <- pvm.mu0(upper, kappa, tol) + pvm.mu0(lower, kappa, tol)
         }
      }     
   }
   return(result)
}

#############################################################
#                                                           #
#   qvonmises function                                      #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: August, 12, 2010                                  #
#   Copyright (C) 2010 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2                                             #
#############################################################

qvonmises <- function(p, mu=circular(0), kappa=NULL, from=NULL, tol = .Machine$double.eps^(0.6), control.circular=list(), ...) {
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
   if (is.null(kappa))
      stop("kappa must be provided")

   zeroPvonmisesRad <- function(x, p, mu, kappa) {
      if (is.na(x)) {    
         y <- NA
      } else {
         y <- integrate(DvonmisesRad, mu=mu, kappa=kappa, lower=0, upper=x)$value - p
      }
      return(y)
   }

   value <- rep(NA, length(p))
   sem <- options()$show.error.messages
   options(show.error.messages=FALSE)
   for (i in 1:length(p)) {
         res <- try(uniroot(zeroPvonmisesRad, p=p[i], mu=mu, kappa=kappa, lower=0, upper=2*pi-epsilon, tol=tol))
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

#############################################################
#                                                           #
#   dmixedvonmises function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 18, 2011                                    #
#   Copyright (C) 2011 Claudio Agostinelli                  #
#                                                           #
#   Version 0.9-2                                           #
#############################################################

dmixedvonmises <- function(x, mu1, mu2, kappa1, kappa2, prop) {
  if (missing(mu1) || length(mu1)!=1)
    stop("the mean direction parameter 'mu1' is mandatory and it must have length 1")
  if (missing(kappa1) || length(kappa1)!=1)
    stop("the concentration parameter 'kappa1' is mandatory and it must have length 1")
  if (missing(mu2) || length(mu2)!=1)
    stop("the mean direction parameter 'mu2' is mandatory and it must have length 1")
  if (missing(kappa2) || length(kappa2)!=1)
    stop("the concentration parameter 'kappa2' is mandatory and it must have length 1")
  if (missing(prop) || length(prop)!=1 || prop > 1 || prop < 0)
    stop("the proportion parameter 'prop' is mandatory and it must have a value between 0 and 1")  
  x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
  mu1 <- conversion.circular(mu1, units="radians", zero=0, rotation="counter")
  mu2 <- conversion.circular(mu2, units="radians", zero=0, rotation="counter")
  mu1 <- as.vector(mu1)
  kappa1 <- as.vector(kappa1)
  mu2 <- as.vector(mu2)
  kappa2 <- as.vector(kappa2)    
  if (kappa1 < 0)
    stop("the concentration parameter 'kappa1' must be non negative")
  if (kappa2 < 0)
    stop("the concentration parameter 'kappa2' must be non negative")
  attr(x, "class") <- attr(x, "circularp") <-  NULL
  attr(mu1, "class") <- attr(mu1, "circularp") <-  NULL
  attr(mu2, "class") <- attr(mu2, "circularp") <-  NULL
    
  DmixedvonmisesRad(x, mu1, mu2, kappa1, kappa2, prop)
}

DmixedvonmisesRad <- function(x, mu1, mu2, kappa1, kappa2, prop) {
  vm <- prop/(2 * pi * besselI(x=kappa1, nu=0, expon.scaled = TRUE)) * (exp(cos(x - mu1) - 1))^kappa1 + (1 - prop)/(2 * pi * besselI(x=kappa2, nu=0, expon.scaled = TRUE)) * (exp(cos(x - mu2) - 1))^kappa2
  return(vm)
}

#############################################################
#                                                           #
#   rmixedvonmises function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 18, 2011                                    #
#   Copyright (C) 2011 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2-6                                           #
#############################################################

rmixedvonmises <- function(n, mu1, mu2, kappa1, kappa2, prop, control.circular=list()) {
  if (missing(mu1) || length(mu1)!=1)
    stop("the mean direction parameter 'mu1' is mandatory and it must have length 1")
  if (missing(kappa1) || length(kappa1)!=1)
    stop("the concentration parameter 'kappa1' is mandatory and it must have length 1")
  if (missing(mu2) || length(mu2)!=1)
    stop("the mean direction parameter 'mu2' is mandatory and it must have length 1")
  if (missing(kappa2) || length(kappa2)!=1)
    stop("the concentration parameter 'kappa2' is mandatory and it must have length 1")
  if (missing(prop) || length(prop)!=1 || prop > 1 || prop < 0)
    stop("the proportion parameter 'prop' is mandatory and it must have a value between 0 and 1")  
  if (is.circular(mu1)) {
    datacircularp <- circularp(mu1)
  } else if  (is.circular(mu2)) {
    datacircularp <- circularp(mu2)
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
   
  mu1 <- conversion.circular(mu1, units="radians", zero=0, rotation="counter")
  mu2 <- conversion.circular(mu2, units="radians", zero=0, rotation="counter")
  mu1 <- as.vector(mu1)
  kappa1 <- as.vector(kappa1)
  mu2 <- as.vector(mu2)
  kappa2 <- as.vector(kappa2)
  if (kappa1 < 0)
    stop("the concentration parameter 'kappa1' must be non negative")
  if (kappa2 < 0)
    stop("the concentration parameter 'kappa2' must be non negative")
  attr(mu1, "class") <- attr(mu1, "circularp") <-  NULL
  attr(mu2, "class") <- attr(mu2, "circularp") <-  NULL

  vm <- RmixedvonmisesRad(n, mu1, mu2, kappa1, kappa2, prop)
  vm <- conversion.circular(circular(vm), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)    
  return(vm)
}

RmixedvonmisesRad <- function(n, mu1, mu2, kappa1, kappa2, prop) {
  result <- rep(NA, n)
  test <- runif(n)
  n1 <- sum(test < prop)
  n2 <- n - n1
  res1 <- RvonmisesRad(n1, mu1, kappa1)
  res2 <- RvonmisesRad(n2, mu2, kappa2)
  result[test < prop] <- res1
  result[test >= prop] <- res2
  return(result)
}

#############################################################
#                                                           #
#   pmixedvonmises function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 18, 2011                                    #
#   Copyright (C) 2011 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################

pmixedvonmises <- function(q, mu1, mu2, kappa1, kappa2, prop, from=NULL, tol = 1e-020) {
  if (missing(mu1) || length(mu1)!=1)
    stop("the mean direction parameter 'mu1' is mandatory and it must have length 1")
  if (missing(kappa1) || length(kappa1)!=1)
    stop("the concentration parameter 'kappa1' is mandatory and it must have length 1")
  if (missing(mu2) || length(mu2)!=1)
    stop("the mean direction parameter 'mu2' is mandatory and it must have length 1")
  if (missing(kappa2) || length(kappa2)!=1)
    stop("the concentration parameter 'kappa2' is mandatory and it must have length 1")
  if (missing(prop) || length(prop)!=1 || prop > 1 || prop < 0)
    stop("the proportion parameter 'prop' is mandatory and it must have a value between 0 and 1")
  q <- conversion.circular(q, units="radians", zero=0, rotation="counter")
  mu1 <- conversion.circular(mu1, units="radians", zero=0, rotation="counter")
  mu2 <- conversion.circular(mu2, units="radians", zero=0, rotation="counter")
  mu1 <- as.vector(mu1)
  kappa1 <- as.vector(kappa1)
  mu2 <- as.vector(mu2)
  kappa2 <- as.vector(kappa2)    
  if (kappa1 < 0)
    stop("the concentration parameter 'kappa1' must be non negative")
  if (kappa2 < 0)
    stop("the concentration parameter 'kappa2' must be non negative")
  attr(q, "class") <- attr(q, "circularp") <-  NULL
  attr(mu1, "class") <- attr(mu1, "circularp") <-  NULL
  attr(mu2, "class") <- attr(mu2, "circularp") <-  NULL
  
  if (is.null(from)) {
    from <- 0
  } else {
    from <- conversion.circular(from, units="radians", zero=0, rotation="counter", modulo="2pi")    
  }  
  mu1 <- (mu1-from)%%(2*pi)
  mu2 <- (mu2-from)%%(2*pi)  
  q <- (q-from)%%(2*pi)
  
  p <- prop*PvonmisesRad(q, mu1, kappa1, tol)+(1-prop)*PvonmisesRad(q, mu2, kappa2, tol)
  return(p)
}

