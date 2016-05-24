#############################################################
#
#   bw.cv.mse.circular function
#   Author: Claudio Agostinelli and Eduardo García Portugués
#   Email: claudio@unive.it
#   date: June, 23, 2011
#   Copyright (C) 2011 Claudio Agostinelli and Eduardo García Portugués
#
#   Version 0.2
#
#############################################################

###References: Hall, Watson and Cabrera, 1984
### Cross validation by MSE ###
bw.cv.mse.circular <- function(x, lower=NULL, upper=NULL, tol = 1e-4, kernel = c("vonmises", "wrappednormal"), K = NULL, min.k = 10) {
  kernel <- match.arg(kernel)
# hmax is removed. Keep this code just for notes
#  if (kernel=="vonmises") {
#    if (is.null(hmax))
#      hmax <- 100
#    if (is.null(lower))
#      lower = 0.01 * hmax
#    if (is.null(upper))
#      upper = hmax
#  } else {
#    if (is.null(hmax))
#      hmax <- 1.144 * sqrt(-2*log(RhoCircularRad(x))) * n^(-1/5) 
#    if (is.null(lower))
#      lower = 0.01 * hmax
#    if (is.null(upper))
#      upper = hmax
#  } 
  if (is.null(upper))
    upper = 50
  if (is.null(lower))
    lower = 0.1
  if ((n <- length(x)) < 2L) 
    stop("need at least 2 data points")  
  x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
  attr(x, "class") <- attr(x, "circularp") <- NULL
  if (!is.numeric(x)) 
    stop("invalid 'x'")
  mse.internal <- function(bw, data) {
    ##bw: bw
    ##data: x
    tone <- integrate(function(z) DensityCircularRad(x=data, z=z, bw=bw, kernel=kernel, K=K, min.k=min.k)^2, lower=0, upper=2*pi, abs.tol=1e-4)$value
    rr <- sapply(1:length(data), function(i) DensityCircularRad(x=data[-i], z=data[i], bw=bw, kernel=kernel, K=K, min.k=min.k))
    ttwo <- 2*sum(rr)/length(data)
    result <- tone-ttwo
    return(result)
  }

  bw <- optimize(function(bw) mse.internal(bw, x), lower=lower, upper=upper, tol=tol, maximum = FALSE)$minimum
  if (bw < lower + tol | bw > upper - tol) 
        warning("minimum occurred at one end of the range")
  return(bw)
}

#############################################################
#
#   bw.cv.ml.circular function
#   Author: Claudio Agostinelli and Eduardo García Portugués
#   Email: claudio@unive.it
#   date: June, 23, 2011
#   Copyright (C) 2011 Claudio Agostinelli and Eduardo García Portugués
#
#   Version 0.2
#
#############################################################

### Cross validation by ML ###
bw.cv.ml.circular <- function(x, lower=NULL, upper=NULL, tol = 1e-4, kernel = c("vonmises", "wrappednormal"), K = NULL, min.k = 10) {
  kernel <- match.arg(kernel)
  if (is.null(upper))
    upper = 50
  if (is.null(lower))
    lower = 0.1
  if ((n <- length(x)) < 2L) 
    stop("need at least 2 data points")  
  x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
  attr(x, "class") <- attr(x, "circularp") <- NULL
  if (!is.numeric(x)) 
    stop("invalid 'x'")
  ml.internal <- function(bw, data) {
    ##bw: bw
    ##data: x
    ss <- sapply(1:length(data), function(i) log(DensityCircularRad(x=data[-i], z=data[i], bw=bw, kernel=kernel, K=K, min.k=min.k)))
    result <- sum(ss)/length(data)
    return(result)
  }

  bw <- optimize(function(bw) ml.internal(bw, x), lower=lower, upper=upper, tol=tol, maximum = TRUE)$maximum
  if (bw < lower + tol | bw > upper - tol) 
        warning("minimum occurred at one end of the range")
  return(bw)
}

#############################################################
#
#   bw.nrd.circular function
#   Author: Claudio Agostinelli and Eduardo García Portugués
#   Email: claudio@unive.it
#   date: July, 22, 2011
#   Copyright (C) 2011 Claudio Agostinelli and Eduardo García Portugués
#
#   Version 0.3
#
#############################################################

###References: Taylor (2008) CSDA formula (7)
bw.nrd.circular <- function(x, lower=NULL, upper=NULL, kappa.est=c("ML","trigmoments"), kappa.bias=FALSE, P=3) {
  if (is.null(upper))
    upper = 50
  if (is.null(lower))
    lower = 0.01
  if ((n <- length(x)) < 2L)
    stop("need at least 2 data points")
  x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
  attr(x, "class") <- attr(x, "circularp") <- NULL
  if (!is.numeric(x))
    stop("invalid 'x'")
  if (is.numeric(kappa.est)) {
    if (length(kappa.est) != 1)
      stop("if 'kappa.est' is numeric, its length must be one")
    kappa <- kappa.est
  } else {
    kappa.est <- match.arg(kappa.est)
    if(kappa.est=="ML"){
      kappa <- MlevonmisesRad(x, mu=NULL, kappa=NULL, bias=kappa.bias)[4]
    } else if(kappa.est=="trigmoments"){
      kappa <- rep(NA, P)
      for (p in 1:P) {
        mup <- TrigonometricMomentRad(x, p, center=FALSE)[1]
        const <- mean(cos(p*x-mup))
        Apzero <- function(x) besselI(x, nu=p, expon.scaled = FALSE)/besselI(x, nu=0, expon.scaled = FALSE) - const
        kappa[p] <- uniroot(f=Apzero, lower=lower, upper=upper)$root
      }
       kappa <- max(kappa)
    } else {
      .NotYetImplemented()
    }
  }
  bw <- (3*n*kappa^2*besselI(x=2*kappa, nu=2, expon.scaled = FALSE)*(4*sqrt(pi)*besselI(x=kappa, nu=0, expon.scaled = FALSE)^2)^(-1))^(2/5)
  return(bw)
}
