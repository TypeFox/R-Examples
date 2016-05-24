#############################################################
#                                                           #
#	wle.vonmises function                               #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 10, 2013                            #
#	Version: 0.4-1                                      #
#                                                           #
#	Copyright (C) 2013 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.vonmises <- function(x, boot=30, group, num.sol=1, raf="HD", smooth, tol=10^(-6), equal=10^(-3), max.iter=500, bias=FALSE, mle.bias=FALSE, max.kappa=500, min.kappa=0.01, use.smooth=TRUE,  alpha=NULL, p=2, verbose=FALSE, control.circular=list()) {

##############################################################
## convoluzione di due vonmises
dvm.convolution <- function(theta, mu1, mu2, kappa1, kappa2) {
  if (kappa1 < 0)
    stop("dvm.convolution: kappa1 must be not negative")
  if (kappa2 < 0)
    stop("dvm.convolution: kappa2 must be not negative")
  kappa <- sqrt(kappa1^2+kappa2^2+2*kappa1*kappa2*cos(theta - (mu1+mu2)))
  ans <- besselI(kappa, 0, expon.scaled=TRUE)/(2*pi*besselI(kappa1, 0, expon.scaled=TRUE)*besselI(kappa2, 0, expon.scaled=TRUE)) * exp(kappa-kappa1-kappa2)
  return(ans)
}

    result <- list()

    if (raf!="HD" & raf!="NED" & raf!="SCHI2") stop("Please, choose the RAF: HD=Hellinger Disparity, NED=Negative Exponential Disparity, SCHI2=Symmetric Chi-squares Disparity")

    if (missing(group)) {
        group <- 0
    }

    if (!is.null(alpha))
        if (alpha==-1)
           p <- Inf
        else 
           p <- (alpha + 1)^(-1)

    # Handling missing values
    x <- na.omit(x)
    size <- length(x)
    if (size==0) {
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

    if (group<2) {
        group <- max(round(size/4),1)
        if (verbose) cat("wle.vonmises: dimension of the subsample set to default value: ",group,"\n")
    }

    if (boot<1) {
        stop("Bootstrap replication not in the range")
    }

    if (!(num.sol>=1)) {
        if (verbose) cat("wle.vonmises: number of solution to report set to 1 \n")
        num.sol <- 1
    }

    if (max.iter<1) {
        if (verbose) cat("wle.vonmises: max number of iteration set to 500 \n")
        max.iter <- 500
    }

    if (tol<=0) {
        if (verbose) cat("wle.vonmises: the accuracy must be positive, using default value: 10^(-6) \n")
        tol <- 10^(-6)
    }

    if (equal<=tol) {
        if (verbose) cat("wle.vonmises: the equal parameter must be greater than tol, using default value: tol+10^(-3)\n")
        equal <- tol+10^(-3)
    }
    if (p< -1) {
        if (verbose) cat("wle.vonmises: the 'p' ('alpha') parameter must be greater than or equal to -1 (0), using default value: p=2 (alpha=-1/2) \n")
        p <- 2
    }

    if (smooth <= 0) {
        if (verbose) stop("wle.vonmises: the smooth parameter must be positive \n")
    }

tot.sol <- 0
not.conv <- 0
iboot <- 0

while (tot.sol < num.sol & iboot < boot) {
   iboot <- iboot + 1
   continue <- TRUE
   start.iter <- 0
   while (continue & start.iter <= max.iter) {
          x.boot <- sample(x, size=group, replace = FALSE)
          temp <- circular:::MlevonmisesRad(x, bias=mle.bias)
          mu <- temp[1]
          kappa <- temp[4]
          continue <- !is.finite(mu) | is.na(mu) | !is.finite(kappa) | is.na(kappa) | kappa > max.kappa | kappa <= 0
          start.iter <- start.iter + 1 
   }
   
   if (start.iter >= max.iter)
       stop("wle.vonmises: the procedure is not able to find a good starting points, please checks the parameters")
     
   if (verbose)
       cat("starting value mu: ", mu, "  kappa: ", kappa, "\n")
   x.diff <- tol + 1
   iter <- 0
   while (x.diff > tol & iter < max.iter) {
       iter <- iter + 1
       mu.old <- mu
       kappa.old <- kappa
### modifica   
       sm <- min(smooth, 1300/kappa)
###       cat("sm ", sm, "\n")
       ff <- circular:::DensityCircularRad(x, z=x, bw=sm*kappa, kernel='vonmises')
       
       if (use.smooth) {
           mm <- dvm.convolution(theta=x, mu1=0, mu2=mu, kappa1=sm, kappa2=kappa)
       } else {
           mm <- circular:::DvonmisesRad(x=x, mu=mu, kappa=kappa)
       }
       if (any(is.nan(mm)) | any(mm==0)) {
           iter = max.iter
       } else { 
           dd <- ff/mm - 1

           if (p==Inf && raf=="HD") {
               ww <- log(dd+1)
           } else {
               ww <- switch(raf,
                 HD =  p*((dd + 1)^(1/p) - 1) ,
	             NED =  2 - (2 + dd)*exp(-dd) ,
	             SCHI2 =  1-(dd^2/(dd^2 +2)) )       
           }
               
           if (raf=="HD" | raf=="NED") {
                ww <- (ww + 1)/(dd + 1)
           }

           ww[ww > 1] <- 1
           ww[ww < 0] <- 0
           sww <- sum(ww)
           wsinr <- ww%*%sin(x)
           wcosr <- ww%*%cos(x)
           mu <- atan2(wsinr, wcosr)
           kappa <- A1inv(ww%*%cos(x - mu)/sww)

           if (bias == TRUE & !is.na(kappa)) {
               if (kappa < 2) {
                   kappa <- max(kappa - 2 * (sww * kappa)^-1, 0)
               } else {
                   kappa <- ((sww - 1)^3 * kappa)/(sww^3 + sww)
               }
           }
           if (!is.na(kappa) && kappa <= 0) {
               if (verbose) cat("wle.vonmises: kappa is set to min.kappa since it is negative \n")             
               kappa <- min.kappa
           }       
           x.diff <- max(abs(mu - mu.old), abs(kappa - kappa.old))
           if (is.na(x.diff)) iter <- max.iter

      }

   }
#### end of while (x.diff > tol & iter < max.iter)

   if (iter < max.iter) {

   if (tot.sol==0) {
      mu.store <- mu
      kappa.store <- kappa
      w.store <- ww
      tot.store <- sum(ww)/size
      f.store <- ff
      m.store <- mm
      d.store <- dd
      tot.sol <- 1
   } else {
      if (min(abs(mu.store%%(2*pi) - mu%%(2*pi))) > equal | min(abs(kappa.store - kappa)) > equal) {
          mu.store <- c(mu.store, mu)
          kappa.store <- c(kappa.store, kappa)
          w.store <- rbind(w.store, ww)
          tot.store <- c(tot.store, sum(ww)/size)
          f.store <- rbind(f.store, ff)
          m.store <- rbind(m.store, mm)
          d.store <- rbind(d.store, dd)
          tot.sol <- tot.sol + 1
      }
   }

   } else not.conv <- not.conv + 1
   
}
##### end of while (tot.sol < num.sol & iboot < boot)

if (tot.sol) {
    result$mu <- conversion.circular(circular(mu.store), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
    result$kappa <- c(kappa.store)   
    result$tot.weights <- tot.store
    result$weights <- w.store
    result$f.density <- f.store
    result$m.density <- m.store
    result$delta <- d.store
    result$tot.sol <- tot.sol
    result$not.conv <- not.conv
} else {
    if (verbose) cat("wle.vonmises: No solutions are fuond, checks the parameters\n")
    result$mu <- conversion.circular(circular(NA), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
    result$kappa <- NA    
    result$tot.weights <- NA
    result$weights <- rep(NA,size)
    result$f.density <- rep(NA,size)
    result$m.density <- rep(NA,size)
    result$delta <- rep(NA,size)
    result$tot.sol <- 0
    result$not.conv <- boot
}

result$call <- match.call()
class(result) <- "wle.vonmises"
return(result)

}

#############################################################
#                                                           #
#	print.wle.vonmises function                         #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 23, 2002                            #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2002 Claudio Agostinelli              #
#                                                           #
#############################################################

print.wle.vonmises <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("mu:\n")
    print.default(format(x$mu, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("kappa:\n")    
    print.default(format(x$kappa, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")    
    cat("\nNumber of solutions ",x$tot.sol,"\n")
    cat("\n")
    invisible(x)
}

