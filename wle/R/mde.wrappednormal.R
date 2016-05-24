#############################################################
#                                                           #
#	mde.wrapperdnormal function                         #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 10, 2013                            #
#	Version: 0.3-2                                      #
#                                                           #
#	Copyright (C) 2013 Claudio Agostinelli              #
#                                                           #
#############################################################

mde.wrappednormal <- function(x, bw, mu=NULL, rho=NULL, sd=NULL, alpha=NULL, p=2, tol=1e-5, n=512, from=circular(0), to=circular(2*pi), lower=NULL, upper=NULL, method="L-BFGS-B", lower.rho=1e-6, upper.rho=1-1e-6, min.sd=1e-3, K=NULL, min.k=10, control.circular=list(), ...) {
   result <- list()
     h.fun <- function(x, xpoints, n, ffty, K, p) {
       k <- c(circular:::DwrappednormalRad(x=xpoints, mu=x[1], rho=x[2], K=K))
       k[k <= 2*.Machine$double.eps] <- 2*.Machine$double.eps
       if (is.finite(p)) {
           k <- k^(1-1/p)
           dist <- Re(2*pi*ffty%*%Conj(fft(k))/(n^2))
           dist <- p^2/(1-p) * (dist - 1)/2
       } else {
           y <- log(k/ffty)
           dist <- Re(2*pi*fft(y)%*%Conj(fft(k))/(n^2))
       }
######       cat(x, dist, "\n")
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
    y <- Re(circular:::DensityCircularRad(x=x, z=z, bw=bw, kernel="wrappednormal", K=K, min.k=min.k))
    y[y <= 2*.Machine$double.eps] <- 2*.Machine$double.eps

    if (is.finite(p)) {
        ffty <- fft(y^(1/p))
    } else {
        ffty <- y
    }
        
    if (is.null(mu) | is.null(rho) | is.null(sd)) {
        res <- circular:::MlewrappednormalRad(x=x, mu=mu, rho=rho, sd=sd, min.sd=min.sd, K=NULL, min.k=min.k, tol=tol, max.iter=100, verbose=FALSE)
        mu <- res[1]
        rho <- res[2]
        sd <- res[3]
    }

    if (rho < 0 | rho > 1)
        stop("rho must be between 0 and 1")

    if (is.null(lower))
       lower <- c(mu - pi, lower.rho)
    if (is.null(upper))
       upper <- c(mu + pi, upper.rho)

    if (is.null(K)) {
        range <- max(mu, x) - min(mu, x)
        K <- (range+6*sd)%/%(2*pi)+1
        K <- max(min.k, K)
    }
     
    res <- try(optim(par=c(mu, rho), fn=h.fun, lower=lower, upper=upper, method=method, xpoints=z, n=n, ffty=ffty, K=K, p=p, ...))
    if (is.list(res) && res$convergence==0) {
        result$dist <- res$value
        result$mu <- res$par[1]
        result$rho <- res$par[2]
        result$sd <- sqrt(-2*log(result$rho))        
        result$k <- circular:::DwrappednormalRad(x=z, mu=result$mu, rho=result$rho, K=K)
    } else {
        result$dist <- NA
        result$mu <- NA
        result$rho <- NA
        result$sd <- NA
        result$k <- rep(NA, length(x))
    }
    result$mu <- conversion.circular(circular(result$mu), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
    result$call <- match.call()
    result$data <- x
    result$x <- z
    result$y <- y
    class(result) <- "mde.wrappednormal"
    return(result)
}

#############################################################
#                                                           #
#	print.mde.wrappednormal function                    #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: June, 11, 2006                                #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2006 Claudio Agostinelli              #
#                                                           #
#############################################################

print.mde.wrappednormal <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("mu:\n")
    print.default(format(x$mu, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("rho:\n")    
    print.default(format(x$rho, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")    
    cat("sd:\n")    
    print.default(format(x$sd, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")    
    invisible(x)
}
