#############################################################
#                                                           #
#	mde.vonmises function                               #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 10, 2013                            #
#	Version: 0.2-3                                      #
#                                                           #
#	Copyright (C) 2013 Claudio Agostinelli              #
#                                                           #
#############################################################

mde.vonmises <- function(x, bw, mu=NULL, kappa=NULL, n = 512, from=circular(0), to=circular(2*pi), lower=NULL, upper=NULL, method="L-BFGS-B", lower.kappa=.Machine$double.eps, upper.kappa=Inf, alpha=NULL, p=2, control.circular=list(), ...) {

  result <- list()
  
  h.fun <- function(x, xpoints, n, ffty, p) {
    k <- c(circular:::DvonmisesRad(x=xpoints, mu=x[1], kappa=x[2]))
    if (is.finite(p)) {
      k <- k^(1-1/p)
      dist <- Re(2*pi*ffty%*%Conj(fft(k))/(n^2))
      dist <- p^2/(1-p) * (dist - 1)/2
    } else {
      y <- log(k/ffty)
      dist <- Re(2*pi*fft(y)%*%Conj(fft(k))/(n^2))
    }
    return(dist)
  }
  if (!is.numeric(from))
    stop("argument 'from' must be numeric")      
  if (!is.numeric(to))
    stop("argument 'to' must be numeric")      
  if (!is.finite(from)) 
    stop("non-finite `from'")
  if (!is.finite(to)) 
    stop("non-finite `to'")
  if (!is.numeric(n))
    stop("argument 'n' must be numeric")
  n <- round(n)
  if (n <=0)
    stop("argument 'n' must be integer and positive")         
  if (!is.numeric(x)) 
    stop("argument 'x' must be numeric")
  if (!is.null(alpha))
    if (alpha==-1)
      p <- Inf
    else 
      p <- (alpha + 1)^(-1)

  if (p < -1) {
    cat("mde.vonmises: the 'p' (alpha) parameter must be greater than or equal to -1 (0), using default value: 2 (-1/2) \n")
    p <- 2
  }

# Handling missing values
  x <- na.omit(x)
  nx <- length(x)
  if (nx==0) {
    warning("No observations (at least after removing missing values)")
    return(NULL)
  }
  if (is.circular(x)) {
    datacircularp <- circularp(x)     
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
  if (!is.null(mu)) {
     mu <- conversion.circular(mu, units="radians", zero=0, rotation="counter")
     attr(mu, "class") <- attr(mu, "circularp") <- NULL
  }
  from <- conversion.circular(from, units="radians", zero=0, rotation="counter")
  attr(from, "class") <- attr(from, "circularp") <- NULL
  to <- conversion.circular(to, units="radians", zero=0, rotation="counter")
  attr(to, "class") <- attr(to, "circularp") <- NULL

  n <- max(n, 512)
  if (n > 512)
    n <- 2^ceiling(log2(n))
  z <- seq(from=from, to=to, length=n)
  y <- Re(circular:::DensityCircularRad(x=x, z=z, bw=bw, kernel="vonmises", K=NULL, min.k=10))
  y[y <= 2*.Machine$double.eps] <- 2*.Machine$double.eps
    
  if (is.finite(p)) {
    ffty <- fft(y^(1/p))
  } else {
    ffty <- y
  }
        
  if (is.null(mu) | is.null(kappa)) {
    res <- circular:::MlevonmisesRad(x=x, mu=mu, kappa=kappa, bias=FALSE)
    mu <- res[1]
    kappa <- res[4]
  }
  if (is.null(lower))
    lower <- c(mu - pi, lower.kappa)
  if (is.null(upper))
    upper <- c(mu + pi, upper.kappa)
  res <- try(optim(par=c(mu, kappa), fn=h.fun, lower=lower, upper=upper, method=method, xpoints=z, n=n, ffty=ffty, p=p, ...))
    
  if (is.list(res) && res$convergence==0) {
    result$dist <- res$value
    result$mu <- res$par[1]
    result$kappa <- res$par[2]
    result$k <- circular:::DvonmisesRad(x=z, mu=result$mu, kappa=result$kappa)
  } else {
    result$dist <- NA
    result$mu <- NA
    result$kappa <- NA
    result$k <- rep(NA, length(z))
  }
  result$mu <- conversion.circular(circular(result$mu), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
  result$call <- match.call()
  result$data <- x
  result$x <- z
  result$y <- y
  class(result) <- "mde.vonmises" 
  return(result)
}

#############################################################
#                                                           #
#	print.mde.vonmises function                         #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: June, 11, 2006                                #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2006 Claudio Agostinelli              #
#                                                           #
#############################################################

print.mde.vonmises <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
  cat("mu:\n")
  print.default(format(x$mu, digits=digits),
		  print.gap = 2, quote = FALSE)
  cat("\n")
  cat("kappa:\n")    
  print.default(format(x$kappa, digits=digits),
		  print.gap = 2, quote = FALSE)
  cat("\n")    
  invisible(x)
}

