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
# but WITHOUT ANY WARRANTY; withx even the implied warranty of
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


## plot.tgp:
##
## generic plot method for tgp-class objects -- handles sensitivity and
## multi-resolution plots through other interfaces after some pre-processing.
## Standard tgp 1-d plots are handled directly, and 2-d projections and slices
## are also handled through other interfaces after a small amount of
## pre-processing

"plot.tgp" <-
function(x, pparts=TRUE, proj=NULL, slice=NULL, map=NULL, as=NULL, 
         center="mean", layout="both", main=NULL, xlab=NULL, ylab=NULL,
         zlab=NULL, pc="pc", gridlen=c(40,40), span=0.1,
         legendloc="topright", maineff=TRUE, mrlayout="both", rankmax=20, ...)
{
  ## check for valid layout
  if(layout != "both" && layout != "surf" && layout != "as" && layout != "sens")
    stop("layout argument must be \"both\", \"surf\", \"as\", or \"sens\"");

  ## check if 'as' plots can be made
  if(x$nn == 0 && (!is.null(as) && (as != "s2" && as != "ks2"))) {
     if(layout == "both") {
	cat("cannot make \"as\" plot since x$nn=0, default to layout = \"surf\"\n")
        layout <- "surf"
     } else if(layout == "as") stop("cannot make \"as\" plot since x$nn=0\n")
  }

  ## sensitivity plots
  if(layout == "sens"){
    if(x$sens$par$ngrid==0){ ## make sure that a sens can be plotted
      cat("Cannot make sensitivity plots without sens.* matrices.\n")
      layout = "both"
    } else { ## plot the sens
      sens.plot(x, legendloc=legendloc, maineff=maineff, ...)
      return(invisible());
    }
  }

  ## plots for multi-resolution tgp
  if(x$params$corr == "mrexpsep"){

    ## the "both" method uses the mr.plot function
    if(mrlayout == "both"){
      mr.plot(x,pparts=pparts, proj=proj, center=center, layout="both",
              main=main, xlab=xlab, ylab=ylab, zlab=zlab, legendloc=legendloc,
              gridlen=gridlen, span=span, ...)
      return(invisible())

      ## whereas the "coarse" and "fine" methods use the regular
      ## tgp plotting methods with some minor changes depending on the res
    } else if(mrlayout == "coarse") {
      xTemp <- x; x <- mr.checkrez(x, res=0)
      if((length(x$Zp.mean)+length(x$ZZ.mean)) < 5)
        stop("Cannot plot 'coarse' with less than 5 predictive locations.\n")
    } else { ## same thing for the fine resolution
      xTemp <- x; x <- mr.checkrez(x, res=1)
      if((length(x$Zp.mean)+length(x$ZZ.mean)) < 5)
        stop("Cannot plot 'fine' with less than 5 predictive locations.\n")
    }
  }

  ## standard tgp plotting
  if(x$d == 1) { # plotting 1d data

    if(layout=="both") par(mfrow=c(1,2), bty="n")
    # else par(mfrow=c(1,1), bty="n")
    
    # construct/get graph labels
    if(is.null(xlab)) xlab <- names(x$X)[1]
    if(is.null(ylab)) ylab <- x$response

    # plot means and errors
    if(layout == "both" || layout == "surf") {

      ## choose mean or median for center
      center <- tgp.choose.center(x, center)
      Z.mean <- center$Z;
      smain <- paste(main, ylab, center$name)
      X <- center$X[,1]
      o <- order(X)

      ## plot the data
      plot(x$X[,1],x$Z, xlab=xlab, ylab=ylab, main=smain,...)

      # plot the center (mean)
      lines(X[o], Z.mean[o], ...)

      ## and 0.5 and 0.95 quantiles
      if(center$name == "kriging mean") {
        Zb.q1 <- Z.mean + 1.96*sqrt(c(x$Zp.ks2, x$ZZ.ks2))
        Zb.q2 <- Z.mean - 1.96*sqrt(c(x$Zp.ks2, x$ZZ.ks2))
      } else {
        Zb.q1 <- c(x$Zp.q1, x$ZZ.q1)
        Zb.q2 <- c(x$Zp.q2, x$ZZ.q2)
      }

      ## add the predictive 90% error-bars
      lines(X[o], Zb.q1[o], col=2, ...)
      lines(X[o], Zb.q2[o], col=2, ...)
      
      # plot parts
      if(pparts & !is.null(x$parts) ) { tgp.plot.parts.1d(x$parts) }
    }
    
    # adaptive sampling plotting
    # first, figure out which stats to plot
    if(layout != "surf") { # && !is.null(as)) {

      ## collect the error statistics that the user has requested
      as <- tgp.choose.as(x, as)
      Z.q <- as$criteria
      X <- as$X

      # then plot them
      o <- order(X[,1]);
      plot(X[o,1], Z.q[o], type="l", ylab=as$name, xlab=xlab,
           main=paste(main, as$name), ...)

      ## plot parts
      if(pparts & !is.null(x$parts)) { tgp.plot.parts.1d(x$parts) }

      ## if improv, then add order too
      if(substr(as$name,1,1) == "I") {
        ranks <- x$improv[,2] <= rankmax
        text(X[ranks,1], Z.q[ranks], labels=x$improv[ranks,2], pos=3, ...)
      }
    }

  } else if(x$d >= 2) { # 2-d plotting
    
    if(x$d == 2 || is.null(slice)) { # 2-d slice projection plot
      tgp.plot.proj(x, pparts=pparts, proj=proj, map=map, as=as, center=center,
                    layout=layout, main=main, xlab=xlab, ylab=ylab, zlab=zlab,
                    pc=pc, gridlen=gridlen, span=span, rankmax=rankmax,
                    ...)
    } else { # 2-d slice plot
      tgp.plot.slice(x, pparts=pparts, slice=slice, map=map, as=as, center=center,
                     layout=layout, main=main, xlab=xlab, ylab=ylab, zlab=zlab,
                     pc=pc, gridlen=gridlen, span=span, ...)    
    }
  } else { ## ERROR
    cat(paste("Sorry: no plot defind for ", x$d, "-d tgp data\n", sep=""))
  }
  ## reset the original tgp object for mr.tgp
  if(x$params$corr == "mrexpsep"){ x <- xTemp }
}
