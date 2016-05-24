bootDR <- function(f, bootmethod="param", niter=1001, silent=TRUE)
{
  if (!inherits(f, "DR"))
    stop("Use only with 'DR' objects")
  if (niter<10) 
    stop("niter must be an integer above 10")
  bootmethod <- match.arg(bootmethod, c("param", "nonparam"))
  
  #simulate bootstrap data
  if (bootmethod == "param") { # parametric bootstrap
    rdistname <- paste("r", f$distname, sep="")
    if (!exists(rdistname, mode="function"))
      stop(paste("The ", rdistname, " function must be defined"))
    rdata <- do.call(rdistname, c(list(n=niter*f$n), as.list(f$estimate)))
    dim(rdata) <- c(f$n, niter)
  }
  else { # non parametric bootstrap
    rdata <- sample(f$data, size=niter*f$n, replace=TRUE)
    dim(rdata) <- c(f$n, niter)
  }
  
  #compute bootstrap estimates
  start <- as.list(f$estimate) #a named vector is no longer is accepted as starting values.
  #if (is.null(f$dots))
  func <- function(iter) {
    res <- try(do.call(fitDR, list(x=rdata[, iter], dist=f$distname, start=start)), silent=silent)
    
    if(class(res)[1] == "try-error")
      return(c(rep(NA, length(start)), 100))
    else
      return(c(res$estimate, res$convergence))
  }
  #else
  #  func <- function(iter) {
  #    res <- do.call(fitDR, c(list(x=rdata[, iter], dist=f$distname, start=start, fix.arg=f$fix.arg), f$dots))
  #    return(c(res$estimate, res$convergence))
  #  }
  owarn <- getOption("warn")
  oerr <- getOption("show.error.messages")
  #print(owarn)
  #print(oerr)
  options(warn=ifelse(silent, -1, 0), show.error.messages=!silent)
  resboot <- sapply(1:niter, func)
  options(warn=owarn, show.error.messages=oerr)
  #print(owarn)
  #print(oerr)
  
  rownames(resboot) <- c(names(f$estimate), "convergence")
  if (length(resboot[, 1])>2) {
    estim <- data.frame(t(resboot)[, -length(resboot[, 1])])
    bootCI <- cbind(apply(resboot[-length(resboot[, 1]), ], 1, median, na.rm=TRUE), 
                    apply(resboot[-length(resboot[, 1]), ], 1, quantile, 0.025, na.rm=TRUE), 
                    apply(resboot[-length(resboot[, 1]), ], 1, quantile, 0.975, na.rm=TRUE))
    colnames(bootCI) <- c("Median", "2.5%", "97.5%")
  }
  else {
    estim <- as.data.frame(t(resboot)[, -length(resboot[, 1])])
    names(estim) <- names(f$estimate)
    bootCI <- c(median(resboot[-length(resboot[, 1]), ], na.rm=TRUE), 
                quantile(resboot[-length(resboot[, 1]), ], 0.025, na.rm=TRUE), 
                quantile(resboot[-length(resboot[, 1]), ], 0.975, na.rm=TRUE)) 
    names(bootCI) <- c("Median", "2.5%", "97.5%") 
  } 
  
  # code of convergence of the optimization function for each iteration
  converg <- t(resboot)[, length(resboot[, 1])]
  
  res <- structure(list(estim=estim, converg=converg, 
                        method=bootmethod, nbboot=niter, CI=bootCI, fitpart=f), 
                   class=c("bootDR", "bootdist"))
  res    
}