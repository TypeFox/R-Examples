#******************************************************************************* 
#
# Bayesian Regression and Adaptive Sampling with Gaussian Process Trees
# Copyright (C) 2005, University of California
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (rbgramacy@ams.ucsc.edu)
#
#*******************************************************************************

## check.sens:
##
## function to check the sens.p argument provided as input
## to the sens or b* functions, depending on the imput dimension
## d

"check.sens" <-
function(sens.p, d)
{
  ## sanity checks
  if(d==1) stop("You don't need sensitivity analysis for a single variable.")
  if(length(sens.p)!=(4*d+3)) stop("bad sens length.")

  ## nn.lhs is 'nm' in the .cc code.
  nn.lhs <- sens.p[1] 
  nn <-nn.lhs*(d+2)
  
  ## The XX matrix is of the correct size for within the .cc code.
  ## This may or may not be necessary.
  ## The first 3 rows contain the LHS parameters to begin with.
  XX <- matrix(rep(0,nn*d),nrow=nn, ncol=d)
  XX[1:2,] <- matrix(sens.p[2:(2*d+1)], nrow=2) ## this is rect

  ## check shape for validity, and copy to XX
  shape <- XX[3,] <- sens.p[(2*d+2):(3*d+1)]
  if(length(shape) != d || !all(shape >= 0)) { 
    print(shape)
    stop(paste("shape should be a non-negative ", d, "-vector", sep=""))
  }

  ## check mode for validity, and copy to XX
  mode <- XX[4,] <- sens.p[(3*d+2):(4*d+1)]
  if(length(mode) != d) {
    print(mode)
    stop(paste("mode should be a ", d, "-vector", sep=""))
  }

  ## check each coordinate of the mode argument
  for(i in 1:d){
    if(mode[i] < XX[1,i] || mode[i] > XX[2,i]){
      stop(paste("mode ", i, " should be within bounds [", 
                 XX[1,i],", ", XX[2,i],"]", sep=""))
    }
  }
  
  ## Create the Main Effect Grid
  ngrid <- sens.p[4*d+2]
  span <- sens.p[4*d+3]
  if((span > 1) || (span < 0)) stop("Bad smoothing span -- must be in (0,1).")
  MEgrid <- matrix(ncol=d, nrow=ngrid)
  for(i in 1:d){ MEgrid[,i] <- seq(XX[1,i], XX[2,i], length=ngrid) }
  
  ## return
  list(nn=nn, nn.lhs=nn.lhs, ngrid=ngrid, span=span, XX=XX,
       MEgrid=MEgrid)
}


## sens:
##
## code for performaing a sensitivity analysis using the specified
## model and nn.lhs LHS re-sampled predictive grid for each of the T
## rounds under a beta prior specified by shape and mode

"sens" <-
function(X, Z, nn.lhs, model=btgp, ngrid=100, span=0.3, BTE=c(3000,8000,10),
         rect=NULL, shape=NULL, mode=NULL, ...)
{
  ## the format for rect is the same as rect in LHS (ncol=2, nrow=d).
  Xnames <- names(X)
  XZ <- check.matrix(X, Z) 
  X <- data.frame(XZ$X);  names(X) <- Xnames;
  Z <- XZ$Z;

  ## process the rect, shape and mode arguments
  d <- ncol(as.matrix(X))
  if(is.null(rect)) rect <- t(apply(as.matrix(X),2,range))
  else if(nrow(rect) != d || ncol(rect) != 2)
    stop(paste("rect should be a ", d, "x2-vector", sep=""))

  ## check the shape LHS parameter vector
  if(is.null(shape)) shape <- rep(1,d)
  else if(length(shape) != d || !all(shape >= 0)) { 
    print(shape)
    stop(paste("shape should be a non-negative ", d, "-vector", sep=""))
  }

  ## check the mode LHS parameter vector
  if(is.null(mode)) mode <- apply(as.matrix(X),2,mean)
  else if(length(mode) != d) {
    print(mode)
    stop(paste("mode should be a ", d, "-vector", sep=""))
  }

  ## check the LHS rectangle in the categorical variable context
  for(i in 1:d){
    if(shape[i]==0){
      if(rect[i,1] != 0 || rect[i,2] != 1){
        print(rect[i,])
        stop(paste("rect must be [0,1] for categorical variables (i=",
                   i,", shape[i]=",shape[i],").", sep=""))
      }
    }
  }
  
  ## build the sens parameter
  sens.p <- c(nn.lhs,t(rect),shape,mode,ngrid,span)

  ## run the b* function (model) with the sens parameter, or otherwise
  ## just return the parameter vector and do nothing
  if(!is.null(model)){ return(model(X,Z,sens.p=sens.p,BTE=BTE,...)) }
  else{ return(sens.p) }  
}


## sens.plot:
##
## function for plotting the results of a sensitivity analysis --
## intended to be used instead of plot.tgp.  The type of plot retulting
## depends on whether main effects are to be plotted or not

"sens.plot" <-
function(s, maineff=TRUE, legendloc="topright",  ...)
{

  ## colors used for each effect (col of X)
  cols = rainbow(s$d)

  ## extract some useful things from the tgp-object 's'
  nom <- names(s$X)
  sens <- s$sens
  Zmean <- sens$ZZ.mean
  Zq1 <- sens$ZZ.q1
  Zq2 <- sens$ZZ.q2

  ## if maineff is logical then the S & T stats will get plotted
  if(is.logical(maineff)){

    ## put X on a mean 0 range 1 scale
    X <- mean0.range1(sens$Xgrid)$X

    ## plot the main effects or not?
    if(maineff){
      par(mfrow=c(1,3), ...)
      X <- mean0.range1(sens$Xgrid)$X

      ## plot each of the main effects in the same window -- start with the 1st
      plot(X[,1], Zmean[,1], main="Main Effects",
           ylab="response", xlab="scaled input",
           col=cols[1], typ="l", lwd=2, ylim=range(as.vector(Zmean)), ...)

      ## and then proceed with the rest
      for(i in 2:s$d){  
        if(nlevels(factor(Zmean[,i]))==3){ ## discrete response ... Taddy is this right?
          segments(-.5, Zmean[1,i], 0, Zmean[1,i], lwd=2, col=cols[i])
          segments(0, Zmean[nrow(Zmean),i], .5, Zmean[nrow(Zmean),i], lwd=2, col=cols[i])
        } else{ lines(X[,i], Zmean[,i], lwd=2, col=cols[i]) } ## continuous response
      }

      ## add a legend to the plot so we can see which colours are for which effects
      legend(x=legendloc, legend = names(s$X), col=cols, fill=cols)
    }
    else{ par(mfrow=c(1,2), ...) }

    ## plot the S and T statistics
    ## S stats first
    boxplot(data.frame(sens$S), names=names(s$X),
            main="1st order Sensitivity Indices", 
	    xlab="input variables", ylab="", ...)

    ## then T stats
    T0 <- sens$T
    T0[sens$T<0] <- 0
    boxplot(data.frame(T0), names=names(s$X),
            main="Total Effect Sensitivity Indices", 
	    xlab="input variables", ylab="", ...)
    
  } else {
    ## only make a main effects plots
    
    ## set up the plot
    X <- sens$Xgrid
    ME <- c(maineff)
    pdim <- dim(as.matrix(maineff))
    par(mfrow=pdim, ...)

    ## for each Main Effect
    for(i in ME){
      ## discrete response ... Taddy is this right?
      if(nlevels(factor(Zmean[,i]))==3){
        plot(c(0,1) ,c(Zmean[1,i],Zmean[nrow(Zmean),i]), main="", ylab="response", 
             xlab=nom[i], col=cols[i], pch=20, cex=2, xlim=c(-.5,1.5), xaxt="n",
             ylim=c(min(Zq1[,i]), max(Zq2[,i])))
        axis(1, at=c(0,1))
        segments(-.1, Zq1[1,i], .1, Zq1[1,i], lwd=2, col=cols[i], lty=2)
        segments(.9, Zq1[nrow(Zq1),i], 1.1, Zq1[nrow(Zq1),i], lwd=2, col=cols[i], lty=2)
        segments(-.1, Zq2[1,i], .1, Zq2[1,i], lwd=2, col=cols[i], lty=2)
        segments(.9, Zq2[nrow(Zq2),i], 1.1, Zq2[nrow(Zq2),i], lwd=2, col=cols[i], lty=2)
      }
      else{ ## continuous response
        plot(X[,i], Zmean[,i], main="", ylab="response", xlab=nom[i],
             col=cols[i], typ="l", lwd=2, 
             ylim=c(min(Zq1[,i]), max(Zq2[,i])), ...)
        lines(X[,i], Zq1[,i], col=cols[i], lty=2)
        lines(X[,i], Zq2[,i], col=cols[i], lty=2)
      }
    }

    ## add a title to the plot
    mtext(text="Main effects: mean and 90 percent interval", line=-2,
          outer=TRUE, font=2)
  }  
}
