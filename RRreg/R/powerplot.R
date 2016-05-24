#' Power plots for multivariate RR methods
#' 
#' Uses the function \code{\link{RRsimu}} to estimate the power of the multivariate RR methods (correlation \code{\link{RRcor}}, logistic regression \code{\link{RRlog}}, and/or linear regression \code{\link{RRlin}}.
#' 
#' @param numRep number of boostrap replications
#' @param n vector of samples sizes
#' @param pi true prevalence
#' @param cor vector of true correlations
#' @param b.log vector of true logistic regression coefficients
#' @param model randomized response model
#' @param p randomization probability
#' @param method multivariate RR method
#' @param complyRates probability of compliance within carriers/noncarriers of sensitive attribute
#' @param sysBias probability of responding 'yes' in case of noncompliance
#' @param groupRatio ratio of subgroups in two-group RR designs
#' @param alpha type-I error used to estimate power
#' @param nCPU number of CPUs to be used
#' @param show.messages toggle printing of progress messages
#' 
#' @return
#' a list of the class \code{powerplot} containing an array \code{res} with the power estimates and details of the simulation (e.g., model, p, pi, etc.)
#' @examples 
#' # Not run
#' # pplot <- powerplot(100, n=c(150,250), cor=c(0,.3,.5),
#' #                   method="RRlog", pi=.6, model="Warner", p=.3)
#' # plot(pplot)
#' @seealso \code{\link{RRsimu}} for Monte-Carlo simulation / parametric bootstrap
#' @export
powerplot <- function(numRep, n=c(100,500,1000), pi, cor=c(0,.1,.3), b.log=NULL, 
                      model, p, method=c("RRcor","RRlog","RRlin"),
                      complyRates = c(1,1), sysBias = c(0,0), 
                      groupRatio=.5, alpha=.05, nCPU=1,show.messages = TRUE){
  # check input n and cor
  if (min(n)<1 || any(floor(n) != n)){
    warning("n must be a vector of positive integers to define the sample sizes")
  }
  if(is.null(b.log)){
    if (min(cor)< -1 || max(cor)>1)
      warning("values for 'cor' must be within the interval [-1,1]")
    ncor <- length(cor)
    b.log <- rep(0, ncor)
  }else{
    warning("Data are generated under the logisitc regression model using the argument 'b.log' (the argument 'cor' is irrelevant).")
    cor <- b.log
    ncor <- length(b.log)
    cor <- rep(0, ncor)
  }
  
  nmet <- length(method)
  res <- array(NA, dim=c(length(n), length(cor), nmet), dimnames=list(n,cor, method))
  
  for (nn in 1:length(n)){
    for (c in 1:ncor){
      sim <- RRsimu(numRep=numRep, n=n[nn], pi=pi, model=model, p=p, method=method,
                    cor=cor[c], b.log=b.log[c], complyRates = complyRates, 
                    sysBias = sysBias, groupRatio=groupRatio, alpha=alpha, 
                    nCPU=nCPU, getPower=T)
      res[nn,c,] <- sim$power
      if(show.messages)
        print(paste("n =",nn, "of", length(n),"; C =",c, "of",length(cor)))
      
      if(! any(is.na(names(sim$power))))
        dimnames(res)[[3]] <- names(sim$power)
    }
    
  }
  

  
  out <- list(res=res, numRep=numRep, n=n, pi=pi, cor=cor, b.log=b.log, 
              model=model, p=p, method=method,
              complyRates = complyRates, sysBias = sysBias, 
              groupRatio=groupRatio, alpha=alpha)
  class(out) <- "powerplot"
  plot(out)
  
  out
}

#' @export
print.powerplot <- function(x, ...){
  cat("Use plot(powerplot(...)) to get a graph.\n")
  print(x$res)
}

#'Plot power of multivariate RR methods
#'
#' Plot estimated power from Monte Carlo simulation as a function of the sample size, separately for different effect sizes and multivariate RR methods
#' @param x a \code{\link{powerplot}} object
#' @param ... ignored
#' 
#' @export
plot.powerplot <- function(x,  ...){
  res <- x$res
  dims <- dim(res)
  dimnam <- dimnames(res)
  n <- x$n
  cor <- x$cor
  b.log <- x$b.log
  method <- x$method
  nmet <- length(method)
  model <- x$model
  p <- x$p
  if(!all(b.log==0)){
    cor <- b.log
    leg <- "b.log ="
  }else{
    leg <- "cor ="
  }
  
  cols <- c("black", "red", "blue","brown","darkgreen","gray","violet","green","darkblue")
  xpd <- par("xpd")
  oma <- par("oma")
  mar <- par("mar")
  
  par(mfrow=c(1,nmet+1), mar=c(4,4,3,1)+.1, oma=c(0,0,0,0))
  for (i in 1:nmet){
    if (length(model) == 1){
      main = paste0( model, "(",paste0(p,collapse=","),")", "-",method[i])
    }else{
      main = paste0( model[1], "(",paste0(p[[1]],collapse=","),")","-",
                    model[2], "(",paste0(p[[2]],collapse=","),")", "-",
                    method[i])
    }
    
    ll <- 1
    plot(n, res[,1,i],type="o", ylim=0:1, ylab=paste0("Power (pi = ",x$pi,")"), 
         xlab="Sample size", bty='L', 
         main=main)
    if (length(cor)>1){
      for (c in 2:length(cor)){
        lines(n, res[,c,i],lty=rep(1:6,4)[c], col=cols[c])
        points(n, res[,c,i],pch=c, col=cols[c] )
        ll[i] <- i
      }  
    }
     abline(h=x$alpha, lty="dotted")
  }
  plot.new()
  par(mar=rep(0,4))
  legend("left", legend=cor, pch=ll, lty=ll, col=cols[1:length(cor)], title=leg)
  par(mfrow=c(1,1) , xpd=xpd, mar=mar, oma=oma)
}
