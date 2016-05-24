# Asymptotically efficient Sobol' indices estimation (Monod 2006, Janon et al. 2012)
#
# Alexandre Janon 2012


sobolEff2 <- function(model = NULL, X1, X2, nboot = 0, conf = 0.95, ...) {
  if ((ncol(X1) != ncol(X2)) | (nrow(X1) != nrow(X2))) {
    stop("The samples X1 and X2 must have the same dimensions")
  }
  p <- ncol(X1)

  X <- X1
  for (i in 1:p) {
    Xb <- X2
    Xb[, i] <- X1[, i]
    X <- rbind(X, Xb) 
  }
  
  x <- list(model = model, X1 = X1, X2 = X2, order = order, nboot = nboot,
            conf = conf, X = X, call = match.call())
  class(x) <- "sobolEff2"
  
  if (! is.null(x$model)) {
    response(x, ...)
    x=tell(x, ...)
  }
  
  return(x)
}


estim.sobolEff2 <- function(data, i=1:nrow(data), estimStd=FALSE, conf=0) {
  d=as.matrix(data[i,],nrow=nrow(data))
  p=ncol(data)
  ind=c()
  meanY=mean(d) #mean(d[,1])
  meanY2=mean(d^2) #mean(d[,1]^2)
  if(estimStd) {
		ind2=c()
		ind2$estimStd=c()
		ind2$CIinf=c()
		ind2$CIsup=c()
		ciScale=qnorm((1+conf)/2)
	}
  for(ii in 1:(ncol(d)-1)) {
	  	  commMean=meanY #.5*(meanY+mean(d[,ii+1]))
		  commSS=meanY2 #.5*(meanY2+mean(d[,ii+1]^2))
		  commMean2=commMean*commMean
		  Vt=commSS-commMean2
        cross=mean(d[,1]*d[,ii+1])
        ind[ii]=(cross-commMean2)/Vt
		  if(estimStd) {
			   U=(d[,1]-commMean)*(d[,ii+1]-commMean)-.5*ind[ii]*((d[,1]-commMean)^2+(d[,ii+1]-commMean)^2)
				ind2$estimStd[ii]=sqrt(var(U)/nrow(data))/Vt
				ind2$CIinf[ii]=ind[ii]-ciScale*ind2$estimStd[ii]
				ind2$CIsup[ii]=ind[ii]+ciScale*ind2$estimStd[ii]
		  }
  }
  if(estimStd) {
		ind2$estim=ind
		return(ind2)
  }
  else {
	  return(ind)
  }
}


tell.sobolEff2 <- function(x, y = NULL, return.var = NULL, ...) {
  p <- ncol(x$X1)
  n <- nrow(x$X1)
  data <- matrix(x$y, nrow = n)
  if (x$nboot == 0) {
    x$S <- estim.sobolEff2(data, 1:n, TRUE, x$conf)
  } else {
    S.boot <- boot(data, estim.sobolEff2, R = x$nboot)
    x$S <- bootstats(S.boot, x$conf, "basic")
  }
  return(x)
}


print.sobolEff2 <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (! is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    if (! is.null(x$S)) {
      cat("\n\n\nSobol indices\n")
		if(!is.null(x$S$estimStd)) {
			df=data.frame(pointEstimate=x$S$estim, stdErr=x$S$estimStd, minCI=x$S$CIinf, maxCI=x$S$CIsup)
			colnames(df)=c("estimate","std. error","min. c.i.","max. c.i.")
			print(df)
		} else {
			print(x$S)
		}
    }
  } else {
    cat("(empty)\n")
  }
}


plot.sobolEff2 <- function(x, ylim = c(0, 1), ...) {
  if (! is.null(x$y)) {
    nodeplot(x$S, ylim = ylim)
  }
}
