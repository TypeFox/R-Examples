#############################################################
#                                                           #
#       wle.ar.ao function                                  #
#	Author: Claudio Agostinelli                             #
#	E-mail: claudio@unive.it                                #
#	Date: December, 20, 2003                            #
#	Version: 0.1-3                                          #
#                                                           #
#	Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

wle.ar.ao <- function(x, x.init, x.seasonal.init, coef, ncoef, ncoef.seasonal, period, sigma2, xreg=NULL, raf, smooth, w.level, verbose=FALSE, ao.list, population.size=20, population.choose=5, elements.random=4, num.max, approx.w=TRUE) {

   nused <- length(x)
   xx <- wle.ar.matrix(x=x, x.init=x.init, x.seasonal.init=x.seasonal.init, ncoef=ncoef, ncoef.seasonal=ncoef.seasonal, period=period, xreg=xreg)
resid <- x - c(xx%*%coef)

   ww <- weights <- .Fortran("wlew",
	as.double(resid), 
	as.integer(nused),
	as.double(resid), 
	as.integer(nused), 
	as.integer(raf),
	as.double(smooth),
	as.double(sigma2),
	totweight=double(1),
	weights=double(nused),
	PACKAGE="wle")$weights

   ao.position <- 0
   pos.temp <- 1:nused
   pos.temp <- pos.temp[order(weights)]
   weights.sort <- sort(weights)
   ao.temp <- weights.sort <= w.level
   pos.temp <- pos.temp[ao.temp]
   ao <- rep(FALSE, nused)
   if (length(pos.temp)) {
       pos.temp <- pos.temp[1:min(length(pos.temp),num.max)]
       ao[pos.temp] <- TRUE
   }
   pos <- (1:nused)[ao]

   if (verbose) {
        cat("We have the following observations under the w.level=",w.level,":\n",pos,"\n")
   }

   if (any(ao)) {
       model.in <- vector(length=0)
       for (i in 1:length(ao.list)) {
            if (all(is.element(ao.list[[i]], pos))) {
            temp <- vector(length=0)
            for (j in 1:length(ao.list[[i]])) {
                 temp <- c(temp,(1:length(pos))[pos==ao.list[[i]][j]])
            }
            model.in <- c(model.in, sum(2^(temp-1)))
       }
   }

   num.model <- max(length(model.in),population.size)
   num.pos <- (2^sum(ao))-1
   dim.dim <- floor(log(num.pos,2))+1
   w.tilde <- rep(0,num.model)

   model.in <- c(model.in, sample(x=(1:num.pos), size=(num.model-length(model.in)), replace=TRUE))

   for (isearch in 1:num.model) {
        pos.ao <- sort(pos[binary(model.in[isearch],dim.dim)$dicotomy])
        num.ao <- length(pos.ao)
        x.ao <- x
        xx.ao <- xx
        for (t in pos.ao) {
             if (ncoef) {
                 x.ao[t] <- xx.ao[t,]%*%coef
                 for (tt in 1:ncoef) {
                      if ((t+tt)<=nused) {  
                          xx.ao[t+tt,tt] <- x.ao[t]
                      }
                 }
             }
             if (ncoef.seasonal) {
                 for (tt in 1:ncoef.seasonal) {
                      if ((t+tt*period)<=nused) {  
                          xx.ao[t+tt*period,tt+ncoef] <- x.ao[t]
                      }
                 }
             }
        }
        resid.ao <- x - c(xx.ao%*%coef)
        resid.ao <- resid.ao[-pos.ao]

        if (approx.w) {
            weights.temp <- approx(x=resid, y=ww, xout=resid.ao)$y 
            weights.temp[is.na(weights.temp)] <- 0
            weights.temp[weights.temp > 1] <- 1
            weights.temp[weights.temp < 0] <- 0            
        } else {
            weights.temp <- wle.weights(x=resid.ao, smooth=smooth, sigma2=sigma2, raf=raf, location=TRUE)$weights
        }
        w.tilde[isearch] <- sum(weights.temp)/nused
   }
   model.in <- model.in[order(w.tilde)]
   w.tilde <- sort(w.tilde)
   while ((model.in[1]-model.in[num.model])!=0) {
          num.model.sel <- population.choose
          cum.wtilde <- cumsum(w.tilde)[num.model.sel:num.model]
          pos.child <- vector(length=0)
          while (length(pos.child)==0) {
                 temp <- runif(2,0,cum.wtilde[length(cum.wtilde)])
                 pos.aaa <- min((num.model.sel:num.model)[cum.wtilde > temp[1]])
                 pos.bbb <- min((num.model.sel:num.model)[cum.wtilde > temp[2]])
                 pos.aa <- pos[binary(model.in[pos.aaa],dim.dim)$dicotomy]
                 pos.bb <- pos[binary(model.in[pos.bbb],dim.dim)$dicotomy]
                 pos.child <- c(pos.aa,pos.bb,pos[sample(x=(1:length(pos)), size=elements.random, replace=TRUE)])
                 pos.child <- pos.child[as.logical(sample(x=c(0,1), size=length(pos.child), replace=TRUE))]
                 pos.child <- sort(unique(pos.child))
          }
          temp.child <- vector(length=0)
          for (i in 1:length(pos.child)) {
               temp.child <- c(temp.child,(1:length(pos))[pos==pos.child[i]])
          }
          model.child <- sum(2^(temp.child-1))
          num.child <- length(pos.child)
          x.ao <- x
          xx.ao <- xx
          for (t in pos.child) {
               x.ao[t] <- xx.ao[t,]%*%coef
               if (ncoef) {
                   for (tt in 1:ncoef) {
                        if ((t+tt)<=nused) {  
                            xx.ao[t+tt,tt] <- x.ao[t]
                        }
                   }
               }
               if (ncoef.seasonal) {
                   for (tt in 1:ncoef.seasonal) {
                        if ((t+tt*period)<=nused) {  
                            xx.ao[t+tt*period,tt+ncoef] <- x.ao[t]
                        }
                   }
               }  
          }
          resid.ao <- x - c(xx.ao%*%coef)
          resid.ao <- resid.ao[-pos.child]
          if (approx.w) {
              weights.temp <- approx(x=resid, y=ww, xout=resid.ao)$y 
              weights.temp[is.na(weights.temp)] <- 0
              weights.temp[weights.temp > 1] <- 1
              weights.temp[weights.temp < 0] <- 0
          } else {
              weights.temp <- wle.weights(x=resid.ao, smooth=smooth, sigma2=sigma2, raf=raf, location=TRUE)$weights
          }
          w.tilde.child <- sum(weights.temp)/nused
          w.tilde <- c(w.tilde,w.tilde.child)
          model.in <- c(model.in,model.child)
          model.in <- model.in[order(w.tilde)][-1]
          w.tilde <- sort(w.tilde)[-1]
   }
   if (max(w.tilde)<(sum(weights)/nused)) {
       ao.position <- NULL
   } else {
       ao.position <- sort(pos[binary(model.in[1],dim.dim)$dicotomy])
   }
} else {
   ao.position <- NULL
}

x.ao <- x
xx.ao <- xx
for (t in ao.position) {
     x.ao[t] <- xx.ao[t,]%*%coef
     if (ncoef) {
         for (tt in 1:ncoef) {
              if ((t+tt)<=nused) {  
                  xx.ao[t+tt,tt] <- x.ao[t]
              }
         }
     }
     if (ncoef.seasonal) {
         for (tt in 1:ncoef.seasonal) {
              if ((t+tt*period)<=nused) {  
                  xx.ao[t+tt*period,tt+ncoef] <- x.ao[t]
              }
         }
     }    
}

resid.ao <- x - c(xx.ao%*%coef)
w.temp <- wle.weights(x=resid.ao, smooth=smooth, sigma2=sigma2, raf=raf, location=TRUE)
resid.ao <- resid.ao - w.temp$location

if (verbose) {
    cat("Additive outliers: \n", ao.position, "\n")
}

return(list(x.ao=x.ao, resid.ao=resid.ao, ao.position=ao.position))
}

#############################################################
#                                                           #
#	wle.ar function                                     #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: January, 5, 2004                              #
#	Version: 0.1-5                                      #
#                                                           #
#	Copyright (C) 2004 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.ar <- function(x, order=c(1,0), seasonal=list(order=c(0,0), period=NA), group, group.start, group.step=group.start, xreg=NULL, include.mean=TRUE, na.action=na.fail, tol=10^(-6), tol.step=tol, equal=10^(-3), equal.step=equal, raf="HD", smooth=0.0031, smooth.ao=smooth,  boot=10, boot.start=10, boot.step=boot.start, num.sol=1, x.init=0, x.seasonal.init=0, max.iter.out=20, max.iter.in=50, max.iter.start=200, max.iter.step=500, verbose=FALSE, w.level=0.4, min.weights=0.5, population.size=10, population.choose=5, elements.random=2, wle.start=FALSE, init.values=NULL, num.max=NULL, num.sol.step=2, approx.w=TRUE) {

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

result <- list()
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

d <- order[2]
d.s <- seasonal$order[2]

#if (d | d.s) {
#    if (any(.packages(all.available=TRUE)=="ts")) {
#        library(ts)
#    } else {
#        stop("For Integrated model you need function diff in package ts")
#    }
#}

if(d) x <- diff(x, 1, d)
if(d.s) x <- diff(x, seasonal$period, d.s)

xtsp <- tsp(x)
tsp(x) <- NULL
nd <- d+d.s
nused <- length(x)

if(is.null(xreg)) {
    ncxreg <- 0
} else {
    if(NROW(xreg) != nused) stop("lengths of x and xreg do not match")
        ncxreg <- NCOL(xreg)
}

class(xreg) <- NULL

if (include.mean && (nd==0)) {
    if (is.matrix(xreg) && is.null(colnames(xreg)))
        colnames(xreg) <- paste("xreg", 1:ncxreg, sep="")
    xreg <- cbind(intercept=rep(1, nused), xreg=xreg)
    ncxreg <- ncxreg + 1
}

if (ncxreg) {
    xreg <- as.matrix(xreg)
    if (qr(xreg)$rank < ncol(xreg)) stop("xreg is collinear")
}

if (ncxreg) {
    if (d) xreg <- diff(xreg, 1, d)
    if (d.s) xreg <- diff(xreg, seasonal$period, d.s)
    xreg <- as.matrix(xreg)
    if (qr(xreg)$rank < ncol(xreg)) stop("xreg is collinear")
}

if (is.null(num.max)) num.max <- nused

# start bootstrap iteration
first.time <- TRUE
iboot <- 1
tot.sol <- 0
not.conv <- 0

while(iboot<=boot & tot.sol<num.sol) {
   pos.iboot <- round(runif(1,(group+1),nused))
   x.boot <- x[(pos.iboot-group+1):pos.iboot]
   if (is.matrix(xreg))
       xreg.boot <- xreg[(pos.iboot-group+1):pos.iboot,]
   else
       xreg.boot <- xreg[(pos.iboot-group+1):pos.iboot]
   
   if (!is.null(init.values)) {
       temp <- list(conv=TRUE)
       temp.n <- length(init.values)
       temp$coef <- init.values[1:(temp.n-1)]
       temp$sigma2 <- init.values[temp.n]
   } else {
       temp <- wle.ar.start(x=x.boot, x.init=x.init, x.seasonal.init=x.seasonal.init, ncoef=ncoef, ncoef.seasonal=ncoef.seasonal, period=period, xreg=xreg.boot, raf=raf, smooth=smooth, group=group.start, boot=boot.start, max.iter=max.iter.start, wle.start=wle.start, verbose=verbose)
   }
    
   if (temp$conv) {
       coef.init <- temp$coef
       sigma2.init <- temp$sigma2

       xx <- wle.ar.matrix(x=x, x.init=x.init, x.seasonal.init=x.seasonal.init, ncoef=ncoef, ncoef.seasonal=ncoef.seasonal, period=period, xreg=xreg)
       resid.init <- x - c(xx%*%coef.init)

      if (verbose) {
	  cat("Initial values from the subsample ",iboot,": \n parameters: ", coef.init,"\n sigma2: ",sigma2.init," \n") 
      }

   weights <- .Fortran("wlew",
    	as.double(resid.init),
    	as.integer(nused),
	as.double(resid.init),
	as.integer(nused),
	as.integer(raf),
	as.double(smooth.ao),
	as.double(sigma2.init),
	totweights=double(1),
	weights=double(nused),
	PACKAGE="wle")$weights

   if (sum(weights)/nused >= min.weights) {

      ao.list <- list(0)
      wres <- wle.ar.ao(x=x, x.init=x.init, x.seasonal.init=x.seasonal.init, coef=coef.init, ncoef=ncoef, ncoef.seasonal=ncoef.seasonal, period=period, sigma2=sigma2.init, xreg=xreg, raf=raf, smooth=smooth.ao, w.level=w.level, verbose=verbose, ao.list=ao.list, population.size=population.size, population.choose=population.choose, elements.random=elements.random, num.max=num.max, approx.w=approx.w)

      resid.ao <- wres$resid.ao
      x.ao <- wres$x.ao
      ao.position <- wres$ao.position
      if (!is.null(ao.position)) {
      ao.list <- c(ao.list,list(ao.position))
      }
      ao.position.old <- c(ao.position,0)
      conv <- TRUE
      iter <- 0

      while (!setequal(ao.position,ao.position.old) & conv) {

    	iter <- iter + 1
    	max.tol <- tol + 1
    	ao.position.old <- ao.position	
        iter.int <- 0
    	while(max.tol>tol & conv) {
            iter.int <- iter.int + 1
	    coef.old <- coef.init

	    res <- wle.ar.step(coef=coef.init, ncoef=ncoef, ncoef.seasonal=ncoef.seasonal, period=period, x=x, xreg=xreg, x.init=x.init, x.seasonal.init=x.seasonal.init, raf=raf, smooth=smooth, sigma2=sigma2.init, num.sol=num.sol.step, ao.position=ao.position, group=group.step, boot=boot.step, max.iter=max.iter.step, verbose=verbose, tol=tol.step, equal=equal.step)

	    coef.init <- res$coef
	    sigma2.init <- res$sigma2
            
	    max.tol <- max(abs(coef.old-coef.init))

	    if (iter.int > max.iter.in) {
                if (verbose) cat("Convergence problem: maximum iteration number reached in the inner loop\n")
                conv <- FALSE
            }

	    if(any(!is.finite(coef.init)) | any(!is.finite(sigma2.init)) | (sum(res$weights)/nused < min.weights)) {
	     	if (verbose) {
                    cat("Convergence problem: some values are not finite, bad starting point or sum of the weights less than min.weights\n")
                    cat("Parameters: ",coef.init,"\n")
                    cat("Sigma2: ",sigma2.init,"\n")
                    cat("Sum of weights/size: ",sum(res$weights)/nused,"\n")
                }
	   	conv <- FALSE
	    }
      }
# end while(max.tol>tol & conv)

      if (conv) {      

	    if(verbose) {
	    	cat("iteration: ",iter," \n parameters: ",coef.init," \n sigma2: ",sigma2.init," \n")
	    }

      wres <- wle.ar.ao(x=x, x.init=x.init, x.seasonal.init=x.seasonal.init, coef=coef.init, ncoef=ncoef, ncoef.seasonal=ncoef.seasonal, period=period, sigma2=sigma2.init, xreg=xreg, raf=raf, smooth=smooth.ao, w.level=w.level, verbose=verbose, ao.list=ao.list, population.size=population.size, population.choose=population.choose, elements.random=elements.random, num.max=num.max, approx.w=approx.w)

      resid.ao <- wres$resid.ao
      x.ao <- wres$x.ao
      ao.position <- wres$ao.position
      if (!is.null(ao.position)) {
      ao.list <- c(ao.list,list(ao.position))
      }

      if ((lao <- length(ao.list))>2) {
         for (i.ao in 1:(lao-2)) {
            if (setequal(ao.list[[lao]],ao.list[[i.ao]])) {
                conv <- FALSE
            }
         } 
      }

      if(iter > max.iter.out) {
      if (verbose) cat("Convergence problem: maximum iteration number reached in the outer loop\n")
      conv <- FALSE
      }

    } 

    }
# end while (!setequal(ao.position,ao.position.old) & conv)

    if (conv) {

    xx <- wle.ar.matrix(x=x, x.init=x.init, x.seasonal.init=x.seasonal.init, ncoef=ncoef, ncoef.seasonal=ncoef.seasonal, period=period, xreg=xreg)
    resid.init <- c(x - c(xx%*%coef.init))  
    resid.init <- ts(resid.init, start=start(x), end=end(x), frequency=frequency(x))
    class(resid.init) <- "ts"

    weights.with.ao <- wle.weights(x=resid.init, smooth=smooth, sigma2=res$sigma2, raf=raf, tol=tol, location=TRUE)$weights
    
    resid <- ts(res$resid, start=start(x), end=end(x), frequency=frequency(x))        
    class(resid) <- "ts" 

    xx.ao <- wle.ar.matrix(x=x.ao, x.init=x.init, x.seasonal.init=x.seasonal.init, ncoef=ncoef, ncoef.seasonal=ncoef.seasonal, period=period, xreg=xreg)
    resid.ao <- c(x - c(xx.ao%*%coef.init))  
    resid.ao <- ts(resid.ao, start=start(x), end=end(x), frequency=frequency(x))

    weights.without.ao <- wle.weights(x=resid.ao, smooth=smooth, sigma2=res$sigma2, raf=raf, tol=tol, location=TRUE)$weights
    
    class(resid.ao) <- "ts" 
    x.ao <- ts(res$x.ao, start=start(x), end=end(x), frequency=frequency(x))
    class(x.ao) <- "ts" 

    if(first.time) {
    coef.final <- res$coef
    sigma2.coef.final <- res$sigma2.coef
    weights.final <- res$weights
    weights.final.with.ao <- weights.with.ao
    weights.final.without.ao <- weights.without.ao   
    sigma2.final <- res$sigma2
    resid.final <- resid
    resid.ao.final <- resid.ao
    x.ao.final <- x.ao
    resid.init.final <- resid.init
    ao.position.final <- list(wres$ao.position)
    first.time <- FALSE	
    tot.sol <- 1
    } else {
    if(min(abs(coef.final-res$coef))>equal) {
	tot.sol <- tot.sol+1
	coef.final <- rbind(coef.final,res$coef)
        sigma2.coef.final <- rbind(sigma2.coef.final,res$sigma2.coef)
	weights.final <- rbind(weights.final,res$weights)
        weights.final.with.ao <- rbind(weights.final.with.ao, weights.with.ao)
        weights.final.without.ao <- rbind(weights.final.without.ao, weights.without.ao)        
	sigma2.final <- c(sigma2.final,res$sigma2)
	resid.final <- rbind(resid.final,resid)
        ao.position.final <- c(ao.position.final,list(wres$ao.position))        
	resid.ao.final <- rbind(resid.ao.final,resid.ao)
        x.ao.final <- rbind(x.ao.final,x.ao)
        resid.init.final <- rbind(resid.init.final,resid.init)
	}
    }
    } else {
    	not.conv <- not.conv+1
    }
    } else {
    	not.conv <- not.conv+1
    }
    } else {
    	not.conv <- not.conv+1
    }    

iboot <- iboot+1
}
# end bootstrap iteration

if(tot.sol==0) {
   
    result$coef <- NULL
    result$sigma2.coef <- NULL
    result$sigma2 <- NULL
    result$arma <- arma
    result$resid <- NULL
    result$resid.with.ao <- NULL
    result$resid.without.ao <- NULL
    result$x.ao <- NULL
    result$call <- match.call()
    result$series <- series
    result$weights <- NULL
    result$weights.with.ao <- NULL
    result$weights.without.ao <- NULL   
    result$tot.sol <- tot.sol
    result$not.conv <- not.conv
    result$ao.position <- NULL

} else {

    nm <- NULL
    if(arma[1] > 0) nm <- c(nm, paste("ar", 1:arma[1], sep=""))
    if(arma[2] > 0) nm <- c(nm, paste("ma", 1:arma[2], sep=""))
    if(arma[3] > 0) nm <- c(nm, paste("sar", 1:arma[3], sep=""))
    if(arma[4] > 0) nm <- c(nm, paste("sma", 1:arma[4], sep=""))
    if(ncxreg > 0)
        if(!is.null(cn <- colnames(xreg))) nm <- c(nm, cn)
        else nm <- c(nm, paste("xreg", 1:ncxreg, sep=""))
    if(tot.sol==1) {
    	names(coef.final) <- nm
        names(sigma2.coef.final) <- nm
    } else {
	colnames(coef.final) <- nm
	rownames(coef.final) <- paste("root", 1:tot.sol, sep=" ")
	colnames(sigma2.coef.final) <- nm
	rownames(sigma2.coef.final) <- paste("root", 1:tot.sol, sep=" ")
	rownames(weights.final) <- paste("root", 1:tot.sol, sep=" ")
	rownames(resid.init.final) <- paste("root", 1:tot.sol, sep=" ")
	rownames(resid.final) <- paste("root", 1:tot.sol, sep=" ")
	rownames(resid.ao.final) <- paste("root", 1:tot.sol, sep=" ")
	rownames(x.ao.final) <- paste("root", 1:tot.sol, sep=" ")
    }	
    names(arma) <- c("ar", "ma", "sar", "sma", "period", "diff", "sdiff")


    result$coef <- coef.final
    result$sigma2.coef <- sigma2.coef.final
    result$sigma2 <- sigma2.final
    result$arma <- arma
    result$resid <- resid.final
    result$resid.without.ao <- resid.ao.final
    result$resid.with.ao <- resid.init.final
    result$x <- x
    result$x.ao <- x.ao.final
    result$call <- match.call()
    result$series <- series
    result$weights <- weights.final
    result$weights.with.ao <- weights.final.with.ao
    result$weights.without.ao <- weights.final.without.ao   
    result$tot.sol <- tot.sol
    result$not.conv <- not.conv
    result$ao.position <- ao.position.final
    result$w.level <- w.level
    }

    class(result) <- "wle.arima"
    
return(result)
}

#############################################################
#                                                           #
#	wle.ar.step function                                #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 30, 2003                            #
#	Version: 0.1-3                                      #
#                                                           #
#	Copyright (C) 2003 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.ar.step <- function(coef, ncoef, ncoef.seasonal, period, x, xreg, x.init, x.seasonal.init, raf, smooth, sigma2, ao.position, group, boot, max.iter, verbose, num.sol=2, tol, equal) {

result <- list()
nused <- length(x)

if (is.null(xreg)) {
    ncxreg <- 0
} else {
    ncxreg <- NCOL(xreg)
}
 
x.ao <- x
xx.ao <- wle.ar.matrix(x=x.ao, x.init=x.init, x.seasonal.init=x.seasonal.init, ncoef=ncoef, ncoef.seasonal=ncoef.seasonal, period=period, xreg=xreg)
for (t in ao.position) {
     x.ao[t] <- xx.ao[t,]%*%coef
     if (ncoef) {
         for (tt in 1:ncoef) {
              if ((t+tt)<=nused) {  
                  xx.ao[t+tt,tt] <- x.ao[t]
              }
         }
     }
     if (ncoef.seasonal) {
         for (tt in 1:ncoef.seasonal) {
              if ((t+tt*period)<=nused) {  
                  xx.ao[t+tt*period,tt+ncoef] <- x.ao[t]
              }
         }
     }
}

#     print(coef)
#     print(cbind(ao.position,x.ao[ao.position]))

if (qr(xx.ao)$rank==NCOL(xx.ao)) {
    temp.wle <- wle.lm(x.ao~xx.ao -1, boot=boot, smooth=smooth, num.sol=num.sol, group=group, max.iter=max.iter, tol=tol, equal=equal)
} else {
    if (verbose) cat("wle.ar.step: the matrix is not full rank\n")
    temp.wle <- list()
    temp.wle$tot.sol==0
}

if (temp.wle$tot.sol!=0) {
  
    if (verbose) {
        cat("Number of solutions: ",temp.wle$tot.sol," found on ",num.sol,"\n")
        cat("Parameters: ",temp.wle$coefficients,"\n")
        cat("Sigma2: ",temp.wle$scale^2,"\n")
    }

if (temp.wle$tot.sol>1) {
    ccc <- c(coef,sigma2)
    dist <- rep(0,temp.wle$tot.sol)
    for (k in 1:temp.wle$tot.sol) {
         dist[k] <- sum((ccc-c(temp.wle$coefficients[k,],temp.wle$scale[k]^2))^2)
    }
    root <- (1:temp.wle$tot.sol)[dist==min(dist)]
    if (verbose) cat("We use root: ",root,"\n")
    ttt <- list()
    ttt$coefficients <- temp.wle$coefficients[root,]
    ttt$residuals <- temp.wle$residuals[root,]
    ttt$scale <- temp.wle$scale[root]
    ttt$weights <- temp.wle$weights[root,]
    temp.wle <- ttt
}

   result$coef <- temp.wle$coefficients
   result$resid <- temp.wle$residuals
   result$sigma2 <- temp.wle$scale^2
   result$weights <- temp.wle$weights
   result$weights[ao.position] <- 0 
   result$sigma2.coef <- diag(result$sigma2*solve(t(xx.ao)%*%diag(result$weights)%*%xx.ao, tol=1e-10))
   result$resid.ao <- temp.wle$residuals
   result$x.ao <- x.ao
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
#	wle.ar.matrix function                              #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: September, 26, 2001                           #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.ar.matrix <- function(x, x.init=0, x.seasonal.init=0, ncoef, ncoef.seasonal, period, xreg=NULL) {

nused <- length(x)

if(is.null(xreg)) {
    ncxreg <- 0
    xreg <- NULL
} else {
    ncxreg <- NCOL(xreg)
}

xx <- vector(length=0)
if (length(x.init)==ncoef) {
    x.temp <- c(x.init,x)
} else {
    x.temp <- c(rep(x.init,ncoef),x)
}

for (i in 1:ncoef) {
     xx <- cbind(xx,x.temp[(ncoef-i+1):(nused+ncoef-i)])
}

if (ncoef.seasonal) {
    if (length(x.seasonal.init)==ncoef.seasonal*period) {
        x.temp <- c(x.seasonal.init,x)
    } else {
        x.temp <- c(rep(x.seasonal.init,ncoef.seasonal*period),x)
    }

    for (i in 1:ncoef.seasonal) {
        xx <- cbind(xx,x.temp[((ncoef.seasonal-i)*period+1):(nused+(ncoef.seasonal-i)*period)])
    }
}

if (ncxreg) xx <- cbind(xx,xreg)

return(xx)
}



#############################################################
#                                                           #
#	wle.ar.start function                               #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: September, 26, 2001                           #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.ar.start <- function(x, x.init, x.seasonal.init, ncoef, ncoef.seasonal, period, xreg, raf, smooth, group, boot, max.iter, wle.start, verbose=FALSE) {

result <- list()
nused <- length(x)

xx <- wle.ar.matrix(x=x, x.init=x.init, x.seasonal.init=x.seasonal.init, ncoef=ncoef, ncoef.seasonal=ncoef.seasonal, period=period, xreg=xreg)

if (qr(xx)$rank==NCOL(xx)) {

    if (wle.start) {
        temp <- wle.lm(x~xx -1, raf=raf, smooth=smooth, group=group, boot=boot, max.iter=max.iter, verbose=verbose, num.sol=1)
        tot.sol <- temp$tot.sol
    } else {
        temp <- lm(x~xx -1)
        tot.sol <- 1
    }

    if (tot.sol) {
        resid <- temp$resid
        result$coef <- temp$coef
        result$resid <- temp$resid
        if (wle.start) {
            result$sigma2 <- temp$scale^2
        } else {
            result$sigma2 <- (summary(temp)$sigma)^2
        }
        result$conv <- TRUE
    } else {
        result$coef <- rep(NA,ncoef)
        result$resid <- rep(NA,nused)
        result$sigma2 <- NA
        result$conv <- FALSE
    }
} else {
    result$coef <- rep(NA,ncoef)
    result$resid <- rep(NA,nused)
    result$sigma2 <- NA
    result$conv <- FALSE
}
     
return(result)
}

