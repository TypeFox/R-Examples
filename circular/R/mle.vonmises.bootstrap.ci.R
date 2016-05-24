
#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   mle.vonmises.bootstrap.ci function                      #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: November, 06, 2013                                #
#   Copyright (C) 2013 Claudio Agostinelli                  #
#                                                           #
#   Version 0.3-3                                           #
#############################################################

mle.vonmises.bootstrap.ci <- function(x, mu=NULL, bias = FALSE, alpha = 0.05, reps = 1000, control.circular=list()) {
  
  # Handling missing values
  x <- na.omit(x)
  if (length(x)==0) {
      warning("No observations (at least after removing missing values)")
      return(NULL)
  }
  if (is.circular(x)) {
     datacircularp <- circularp(x)     
  } else
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

  x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
  attr(x, "class") <- attr(x, "circularp") <- NULL
  if (is.null(mu)) {
     sinr <- sum(sin(x))
     cosr <- sum(cos(x))
     mu <- atan2(sinr, cosr)
  } else {
     mu <- conversion.circular(mu, units="radians", zero=0, rotation="counter", modulo="2pi")
     attr(mu, "class") <- attr(mu, "circularp") <- NULL
  }

  result <- MleVonmisesBootstrapCiRad(x, mu, bias, alpha, reps)
  
  result$mu <- conversion.circular(circular(result$mu), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
  result$mu.ci <- conversion.circular(circular(result$mu.ci), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)

  result$call <- match.call()
  result$alpha <- alpha
  class(result) <- "mle.vonmises.bootstrap.ci"
  return(result)
}

MleVonmisesBootstrapCiRad <- function(x, mu, bias, alpha, reps) {
  mean.bs <- boot(data = x, statistic = MleVonmisesMuRad, R = reps, stype="i")
  mean.reps <- mean.bs$t
  mean.reps <- sort(mean.reps %% (2 * pi))
  spacings <- c(diff(mean.reps), mean.reps[1] - mean.reps[reps] + 2 * pi)
  max.spacing <- (1:reps)[spacings == max(spacings)]
  off.set <- 2 * pi - mean.reps[max.spacing + 1]
  if (max.spacing != reps)
     mean.reps2 <- mean.reps + off.set
  else
    mean.reps2 <- mean.reps
  mean.reps2 <- sort(mean.reps2 %% (2 * pi))
  mean.ci <- quantile(mean.reps2, c(alpha/2, 1 - alpha/2))
  if (max.spacing != reps)
    mean.ci <- mean.ci - off.set
      
  kappa.bs <- boot(data = x, statistic = MleVonmisesKappaRad, R = reps, stype="i", mu=mu, bias = bias)
  kappa.reps <- kappa.bs$t
  kappa.ci <- quantile(kappa.reps, c(alpha/2, 1 - alpha/2))
  result <- list()
  result$mu.ci <- mean.ci
  result$mu <- c(mean.reps)
  result$kappa.ci <- kappa.ci
  result$kappa <- c(kappa.reps)
  return(result)
}

MleVonmisesMuRad <- function(x, i) {
   sinr <- sum(sin(x[i]))
   cosr <- sum(cos(x[i]))
   mu <- atan2(sinr, cosr)
   return(mu)
}

MleVonmisesKappaRad <- function(x, i, mu, bias) {
   n <- length(x[i])
   V <- mean(cos(x[i] - mu))
   if (V > 0) {
      kappa <- A1inv(V)
   } else {
      kappa <- 0
   }
   if (bias == TRUE) {
      if (kappa < 2) {
         kappa <- max(kappa - 2 * (n * kappa)^-1, 0)
      } else {
         kappa <- ((n - 1)^3 * kappa)/(n^3 + n)
      }
   }
   return(kappa)
}
      
#############################################################
#                                                           #
#   print.mle.vonmises.bootstrap.ci function                #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: September, 17, 2003                               #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################

print.mle.vonmises.bootstrap.ci <- function(x, ...) {
    cat("Bootstrap Confidence Intervals for Mean Direction and Concentration", "\n")
    cat("Confidence Level:  ", round(100 * (1 - x$alpha),2), "%", "\n")
    cat("Mean Direction:           ", "Low =", round(x$mu.ci[1], 2), "  High =", round(x$mu.ci[2], 2), "\n")
    cat("Concentration Parameter:  ", "Low =", round(x$kappa.ci[1], 2), "  High =", round(x$kappa.ci[2], 2), "\n")

}
        
