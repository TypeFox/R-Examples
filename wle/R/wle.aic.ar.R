#############################################################
#                                                           #
#	wle.aic.ar function                                 #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 24, 2003                            #
#	Version: 0.1-1                                      #
#                                                           #
#	Copyright (C) 2003 Claudio Agostinelli              #
#                                                           #
#############################################################

# La funzione con il metodo WLS potrebbe non funzionare in presenza
# dell'operatore differenza  
# E' una cosa da verificare
 
wle.aic.ar <- function(x, order=c(1,0), seasonal=list(order=c(0,0), period=NA), group, group.start, group.step=group.start, xreg=NULL, include.mean=TRUE, na.action=na.fail, tol=10^(-6), tol.step=tol, equal=10^(-3), equal.step=equal, raf="HD",  var.full=0, smooth=0.0031, smooth.ao=smooth,  boot=10, boot.start=10, boot.step=boot.start, num.sol=1, x.init=0, x.seasonal.init=0, max.iter.out=20, max.iter.in=50, max.iter.start=200, max.iter.step=500, verbose=FALSE, w.level=0.4, min.weights=0.5, population.size=10, population.choose=5, elements.random=2, wle.start=FALSE, init.values=NULL, num.max=NULL, num.sol.step=2, min.weights.aic=0.5, approx.w=TRUE, ask=TRUE, alpha=2, method="WLS") {

if (length(order)!=2) stop("order must have two components")
if (length(seasonal$order)!=2) stop("seasonal$order must have two components")

ncoef <- order[1]
ncoef.seasonal <- seasonal$order[1]

raf <- switch(raf,
	HD = 1,
	NED = 2,
	SCHI2 = 3,
	-1)

if (raf==-1) stop("Please, choose the RAF: HD=Hellinger Disparity, NED=Negative Exponential Disparity, SCHI2=Symmetric Chi-squares Disparity")

series <- deparse(substitute(x))
if(NCOL(x) > 1)
    stop("only implemented for univariate time series")

x <- na.action(as.ts(x))
n <- length(x)
if(is.null(seasonal$period) || is.na(seasonal$period)
    || seasonal$period == 0) seasonal$period <- frequency(x)

period <- seasonal$period
arma <- c(ncoef,0,ncoef.seasonal,0,period,order[2],seasonal$order[2])

if (!(length(x.init)==1 | length(x.init)==ncoef)) stop("x.init must have one or order[1] elements\n")
if (!(length(x.seasonal.init)==1 | length(x.seasonal.init)==ncoef.seasonal*period)) stop("x.seasonal.init must have one or seasonal$order[1]*period elements\n")

pos <- 0
res <- list()

if (include.mean) {
    p.mean <- 1
} else {
    p.mean <- 0
}

for (mm in p.mean:0) {
if (mm) {
    nn <- 0
} else {
    nn <- 1
}
for (ar in ncoef:nn) {
     for (sar in ncoef.seasonal:0) {
          pos <- pos + 1
          pos.order <- c(ar,order[2])
          pos.mean <- mm==1
          pos.seasonal <- list(order=c(sar,seasonal$order[2]), period=seasonal$period)

          if (pos.order[1]==0 & pos.seasonal$order[1]==0 & mm==1) {
              if (method=="WLE") {
                  result <- wle.normal(x=x, group=2, tol=tol, equal=equal, raf=raf, smooth=smooth, boot=boot.start, num.sol=num.sol,  max.iter=max.iter.step, verbose=verbose)
                  result$resid.without.ao <- result$residuals
              } else {
                  result <- list()
                  result$location <- sum(x.ao*weights)/sum(weights)
                  result$scale <- sum((x.ao-result$location)^2*weights)/sum(weights)
                  result$residuals <- x.ao - result$location
                  result$resid.without.ao <- result$residuals
              }    
          } else {
              if (pos==1 | method=="WLE") {
                  result <- wle.ar(x=x, order=pos.order, seasonal=pos.seasonal, group=group, group.start=group.start, group.step=group.step, xreg=xreg, include.mean=pos.mean, na.action=na.fail, tol=tol, tol.step=tol.step, equal=equal, equal.step=equal.step, raf=raf, smooth=smooth, smooth.ao=smooth.ao,  boot=boot, boot.start=boot.start, boot.step=boot.step, num.sol=num.sol, x.init=x.init, x.seasonal.init=x.seasonal.init, max.iter.out=max.iter.out, max.iter.in=max.iter.in, max.iter.start=max.iter.start, max.iter.step=max.iter.step, verbose=verbose, w.level=w.level, min.weights=min.weights, population.size=population.size, population.choose=population.choose, elements.random=elements.random, wle.start=wle.start, init.values=init.values, num.max=num.max, num.sol.step=num.sol.step, approx.w=approx.w)
              } else {
                  result <- wle.ar.wls(ncoef=pos.order[1], ncoef.seasonal=pos.seasonal$order[1], period=pos.seasonal$period, x=x.ao, xreg=xreg, x.init=x.init, x.seasonal.init=x.seasonal.init, weights=weights, verbose=verbose)
              }
          }
              
          res[[pos]] <- list(
                 order=pos.order,
                 seasonal=pos.seasonal,
                 include.mean=pos.mean,           
                 result=result
          )

          if (pos==1) {
              if (res[[1]]$result$tot.sol==0) {
                  stop ("No solutions are found in the full model, please check the parameters")
              }
              if (res[[1]]$result$tot.sol>1) {
                  sum.w <- apply(res[[1]]$result$weights.without.ao,1,sum)
                  num <- apply(res[[1]]$result$x.ao,1,length)
                  good <- min.weights.aic <= sum.w/num
                  if (!any(good)) {
                      stop (paste("Please check min.weight.aic parameter. The sum of weights with respect to the sample size is ",sum.w/num))
                  }
                  if (ask) {
                      tot.sol <- res[[1]]$result$tot.sol
                      temp <- res[[1]]$result
                      cat("Coefficients:\n")
                      print(temp$coeff)
                      cat("Residual variances:\n")                  
                      print(temp$sigma2)
                      cat("Insert the root (a number between 1 and ",tot.sol,")\n")
                      pos.sol <- readline("\n")
                   } else {
                      pos.sol <- min(order(temp$sigma2)[good])
                   }
                   sigma2 <-  res[[1]]$result$sigma2[pos.sol]
                   weights <- res[[1]]$result$weights.without.ao[pos.sol,]
                   x.ao <- res[[1]]$result$x.ao[pos.sol,]
              } else {
                   sigma2 <-  res[[1]]$result$sigma2
                   weights <- res[[1]]$result$weights.without.ao
                   x.ao <- res[[1]]$result$x.ao
              }
          }
     }
}
}

num.model <- length(res)

if (var.full==0) var.full <- sigma2

tot.w <- sum(weights)
waic <- vector(length=0)
fix <- log(2*pi*var.full)*tot.w
for (imodel in 1:num.model) {
     isol <- res[[1]]$result$tot.sol
     ar <- res[[imodel]]$order[1]
     sar <- res[[imodel]]$seasonal$order[1]
     mm <- as.numeric(res[[imodel]]$include.mean)

     if (isol==1) { 
         temp <- sum(weights*res[[imodel]]$result$resid.without.ao^2)/var.full + alpha*(ar+sar+mm) + fix
         waic <- rbind(waic,c(ar,sar,mm,temp))
     } else {
         for (iter in 1:isol) {
              temp <- sum(weights*(res[[imodel]]$result$resid.without.ao[iter,])^2)/var.full + alpha*(ar+sar+mm) + fix
              waic <- rbind(waic,c(ar,sar,mm,temp))
         }
     }
}

if (num.model>1) { 
    colnames(waic) <- c("ar", "sar", "include.mean", "waic")
} else {
    names(waic) <- c("ar", "sar", "include.mean", "waic")
}

result <- list(full.model=res[[1]]$result, waic=waic)          
result$call <- match.call()


class(result) <- "wle.aic.ar"

return(result)         
}


wle.ar.wls <- function(ncoef, ncoef.seasonal, period, x, xreg, x.init, x.seasonal.init, ao.position, weights, verbose) {
  
result <- list()
nused <- length(x)

if (is.null(xreg)) {
    ncxreg <- 0
} else {
    ncxreg <- NCOL(xreg)
}
 
xx.ao <- wle.ar.matrix(x=x, x.init=x.init, x.seasonal.init=x.seasonal.init, ncoef=ncoef, ncoef.seasonal=ncoef.seasonal, period=period, xreg=xreg)

if (qr(xx.ao)$rank==NCOL(xx.ao)) {
    temp.wle <- lm(x~xx.ao -1, weights=weights)
    temp.wle$tot.sol <- 1
  } else {
    if (verbose) cat("wle.ar.wls: the matrix is not full rank\n")
    temp.wle <- list()
    temp.wle$tot.sol <- 0    
}

if (temp.wle$tot.sol!=0) {
  
    if (verbose) {
        cat("Number of solutions: 1\n")
        cat("Parameters: ",temp.wle$coefficients,"\n")
        cat("Sigma2: ",summary(temp.wle)$sigma^2,"\n")
    }
    
   result$coef <- temp.wle$coefficients
   result$resid <- temp.wle$residuals
   result$sigma2 <- summary(temp.wle)$sigma^2
   result$weights <- weights
   result$sigma2.coef <- diag(result$sigma2*solve(t(xx.ao)%*%diag(result$weights)%*%xx.ao, tol=1e-10))
   result$resid.without.ao <- temp.wle$residuals
   result$x.ao <- x
   result$sol <- 1
} else {
   result$sol <- 0
   result$coef <- rep(NA,ncol(xx.ao))
   result$sigma2 <- NA
}

return(result)
}


#############################################################
#                                                           #
#	summary.wle.aic.ar function                         #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: April, 10, 2005                               #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2005 Claudio Agostinelli              #
#                                                           #
#############################################################

summary.wle.aic.ar <- function (object, num.max=20, verbose=FALSE, ...) {

if (num.max<1) {
    if (verbose) cat("summary.wle.aic.ar: num.max can not less than 1, num.max=1 \n")
    num.max <- 1
}

ans <- list()
waic <- object$waic
if (is.null(nmodel <- nrow(waic))) nmodel <- 1
num.max <- min(nmodel,num.max)
if (nmodel!=1) { 
    nvar <- ncol(waic)-1
    waic <- waic[order(waic[,(nvar+1)]),]
    waic <- waic[1:num.max,]
}

ans$waic <- waic
ans$num.max <- num.max
ans$call <- object$call

class(ans) <- "summary.wle.aic.ar"
return(ans)
}

#############################################################
#                                                           #
#	print.wle.aic.ar function                           #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: April, 10, 2005                               #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2005 Claudio Agostinelli              #
#                                                           #
#############################################################

print.wle.aic.ar <- function (x, digits = max(3, getOption("digits") - 3), num.max=max(1, nrow(x$waic)), ...) {
    res <- summary.wle.aic.ar(object=x, num.max=num.max, ...)
    print.summary.wle.aic.ar(res, digits=digits, ...)
}

#############################################################
#                                                           #
#	print.summary.wle.aic.ar function                   #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: April, 10, 2005                               #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2005 Claudio Agostinelli              #
#                                                           #
#############################################################

print.summary.wle.aic.ar <- function (x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")

    cat("\nWeighted Akaike Information Criterion (WAIC):\n")
    if(x$num.max>1) {
    nvar <- ncol(x$waic)-1
    x$waic[,(nvar+1)] <- signif(x$waic[,(nvar+1)],digits)
    } else {
    nvar <- length(x$waic)-1
    x$waic[(nvar+1)] <- signif(x$waic[(nvar+1)],digits)
    }
    print(x$waic)
    cat("\n")

    cat("Printed the first ",x$num.max," best models \n") 
    invisible(x)
}















