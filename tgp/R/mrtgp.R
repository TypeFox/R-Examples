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

## mr.plot:
##
## plotting function for multiresolution tgp-class objects
## (i.e., those with corr=="mrexpsep") -- called by plot.tgp

"mr.plot" <-
function(x, pparts=TRUE, proj=NULL, center="mean", layout="both",
         main=NULL, xlab=NULL, ylab=NULL, zlab=NULL, legendloc="topright",
         gridlen=c(40,40), span=0.1, ...)
{
  ## 1-d plot of 1-d data described by two columns (resolutions)
  if( x$d==2 ){

    ## create plot window
    par(mfrow=c(1,1))

    ## construct axis x&y labels
    if(is.null(xlab)){xlab <- names(x$X)[2]}
    if(is.null(ylab)){ylab <- x$response}

    ## collect the input and predictive data and pred outputs
    center <- tgp.choose.center(x, center)
    o <- order(center$X[,2])
    X <- center$X[o,]
    Z <- center$Z[o]
    smain <- paste(main, ylab, center$name)

    ## collect quantiles
    Z.q1 <- c(x$Zp.q1, x$ZZ.q1)[o]
    Z.q2 <- c(x$Zp.q2, x$ZZ.q2)[o]

    ## plot the coarse and fine input data
    plot(x$X[x$X[,1]==0,2],x$Z[x$X[,1]==0], ylim=range(c(Z,x$Z)),
         xlab=xlab, ylab=ylab, main=smain, col=4)
    lines(x$X[x$X[,1]==1,2],x$Z[x$X[,1]==1], type="p", pch=20,
          col=2)

    ## add a legend
    if(! is.null(legendloc)) 
      legend(legendloc, lty=c(1,2,1,2), col=c("blue", "blue", "red", "red"),
             c(paste("coarse", center$name), "coarse 90% CI",
               paste("fine", center$name), "fine 90% CI"))

    ## extract the coarse and fine resolutions
    f<-X[,1]==1
    c<-X[,1]==0

    ## add the coarse and fine mean and quantiles
    lines(X[c,2], Z[c], col=4)
    lines(X[f,2], Z[f], col=2)
    lines(X[f,2], Z.q1[f], col=2, lty=3)
    lines(X[f,2], Z.q2[f], col=2, lty=3)
    lines(X[c,2], Z.q1[c], col=4, lty=3)
    lines(X[c,2], Z.q2[c], col=4, lty=3)
    if(pparts) tgp.plot.parts.1d(x$parts[,2])
   
   } else { ## make a projection for data is >= 2-d
 
     ## create plot window
     par(mfrow=c(1,2))
     if(is.null(proj)) proj <- c(1,2)

     ## create axis lables -- augment proj argument by one column
     proj <- proj+1
     
     if(is.null(xlab)){xlab <- names(x$X)[proj[1]]}
     if(is.null(ylab)){ylab <- names(x$X)[proj[2]]}

     ## collect the input and predictive data and pred outputs
     ## this plot only plots the mean or median, no errors
     center <- tgp.choose.center(x, center)
     X <- center$X; Z <- center$Z

     ## separate X and Z into coarse and fine
     c<-X[,1]==0; f<-X[,1]==1
     Xc <- as.data.frame(X[c,proj])
     Xf <- as.data.frame(X[f,proj])
     Zc <- Z[c]; Zf <- Z[f]

     ## initialize the projection vectors p*
     nXc <- nrow(Xc);  pc <- seq(1,nXc)
     nXf <- nrow(Xf);  pf <- seq(1,nXf)
     dX <- nrow(X)

     ## plot the coarse predictive (mean or median) surface
     smain <- paste(main, x$response, "coarse", center$name)
     slice.image(Xc[,1], Xc[,2], p=pc, z=Zc, xlab=xlab, ylab=ylab,
                 main=smain, gridlen=gridlen,span=span,
                 xlim=range(X[,proj[1]]), ylim=range(X[,proj[2]]), ...)
     ## add inputs and predictive locations
     points(x$X[x$X[,1]==0,proj], pch=20, ...)
     points(x$XX[x$XX[,1]==0,proj], pch=21, ...)
     
     # plot parts
     if(pparts & !is.null(x$parts)) { tgp.plot.parts.2d(x$parts, dx=proj)}

     ## plot the fine predictive (mean or median) surface
     smain <- paste(main, x$response, "fine", center$name)
     slice.image(Xf[,1], Xf[,2], p=pf, z=Zf, xlab=xlab, ylab=ylab,
                 main=smain, gridlen=gridlen, span=span,
                 xlim=range(X[,proj[1]]), ylim=range(X[,proj[2]]), ...)
     ## add inputs and predictive locations
     points(x$X[x$X[,1]==1,proj], pch=20, ...)
     points(x$XX[x$XX[,1]==1,proj],pch=21, ...)
     
     # plot parts
     if(pparts & !is.null(x$parts)) { tgp.plot.parts.2d(x$parts, dx=proj)}
   }
}


## mr.checkrez:
##
## used for extreacting the predictive surface information for
## one of the two resolutions so that the surface for that
## resolution can be plotted using the regualr tgp plotting
## machinery in plot.tgp

"mr.checkrez" <-
function(b, res)
{
  ## select input data at the desired resolution
  b$d <- b$d-1
  rdata <- b$X[,1]==res
  b$n <- sum(rdata)
  cnames=names(b$X)[-1]
  b$X <- as.data.frame(b$X[rdata,-1])
  colnames(b$X) <- cnames
  b$Z <- b$Z[rdata]

  ## predictive data at input locations for the desired resolution
  b$Zp.mean <- b$Zp.mean[rdata]
  b$Zp.km <- b$Zp.km[rdata]
  b$Zp.q <- b$Zp.q[rdata]
  b$Zp.s2 <- b$Zp.s2[rdata]
  b$Zp.ks2 <- b$Zp.ks2[rdata]
  b$Zp.q1 <- b$Zp.q1[rdata]
  b$Zp.q2<- b$Zp.q2[rdata]
  b$Zp.med <- b$Zp.med[rdata]

  ## predictive data at the predictive locations for the desired resolution
  rpred <- b$XX[,1]==res
  b$nn <- sum(rpred)
  b$XX <- as.data.frame(b$XX[rpred,-1])
  colnames(b$XX) <- cnames
  b$ZZ <- b$ZZ[rpred]
  b$ZZ.mean <- b$ZZ.mean[rpred]
  b$ZZ.km <- b$ZZ.km[rpred]
  b$ZZ.q <- b$ZZ.q[rpred]
  b$ZZ.s2 <- b$ZZ.s2[rpred]
  b$ZZ.ks2 <- b$ZZ.ks2[rpred]
  b$ZZ.q1 <- b$ZZ.q1[rpred]
  b$ZZ.q2<- b$ZZ.q2[rpred]
  b$ZZ.med <- b$ZZ.med[rpred]
  b$improv <- b$improv[rpred,]

  b$parts <- b$parts[,-1]
  
  return(b)
}
