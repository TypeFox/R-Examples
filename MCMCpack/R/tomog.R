##########################################################################
## Tomography Plots for Ecological Inference 
##
## produces tomography plots (see King, 1997, A Solution to the
## Ecological Inference Problem, Princeton University Press)
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## KQ 11/9/2002
##
## Modification added suggested by David Hugh-Jones 6/10/2006
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################


"tomogplot" <-
  function(r0, r1, c0, c1, 
           xlab="fraction of r0 in c0 (p0)",
           ylab="fraction of r1 in c0 (p1)",
           bgcol="white", ...) {
    if (length(r0) != length(r1)) {
      stop("r0 and r1 different lengths in tomogplot().\n")
    }
    if (length(r0) != length(c0)) {
      stop("r0 and c0 different lengths in tomogplot().\n")
    }
    if (length(r0) != length(c1)) {
      stop("r0 and c1 different lengths in tomogplot().\n")
    }
    
    intercept <-  c0/r1
    slope <- -1 * r0/r1
    N <- length(r0)
    
    par(pty="s")
    plot(0:1, 0:1, type="n", main="", xlab=xlab, ylab=ylab)
    rect(0, 0, 1, 1, col=bgcol, lty=0)
    
    for (year in 1:N) {
       if (is.finite(intercept[year]) & is.finite(slope[year]))
			abline(intercept[year], slope[year])
		else abline(v=c0[year]/(c0[year]+c1[year]))		       
    }
    
    rect(-0.05, -0.05, 1.05, 0, col="white", lty=0)
    rect(-0.05, -0.05, 0, 1.05, col="white", lty=0)
    rect(-0.05, 1, 1.05, 1.05, col="white", lty=0)
    rect(1, -0.05, 1.05, 1.05, col="white", lty=0)
    box()
    par(pty="m")
    return(0)
  }

## produces temporally organized tomography plots
## (see King, 1997, A Solution to the Ecological Inference
## Problem, Princeton University Press)
##
## KQ 11/9/2002
## Modification added suggested by David Hugh-Jones 6/10/2006

"dtomogplot" <-
  function(r0, r1, c0, c1, time.vec=NA, delay=0, 
           xlab="fraction of r0 in c0 (p0)",
           ylab="fraction of r1 in c0 (p1)",
           color.palette=heat.colors,
           bgcol="black", ...) {
    if (length(r0) != length(r1)){
      stop("r0 and r1 different lengths in dtomogplot().\n")
    }
    if (length(r0) != length(c0)){
      stop("r0 and c0 different lengths in dtomogplot().\n")
    }
    if (length(r0) != length(c1)){
      stop("r0 and c1 different lengths in dtomogplot().\n")
    }
    if (length(r0) != length(time.vec) & !is.na(time.vec)[1]){
      stop("r0 and time.vec different lengths in dtomogplot().\n")
    }
  
    intercept <-  c0/r1
    slope     <- -1 * r0/r1
    N <- length(r0)
    if (is.na(time.vec)[1])
      time.vec <- 1:N
    col.vec <- color.palette(N)

    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    on.exit(par(par.orig))
    w <- (3 + mar.orig[2]) * par("csi") * 2.54
    layout(matrix(c(2,1), ncol=2), widths=c(1,lcm(w)))
    par(las=1)
    mar <- mar.orig
    mar[4] <- mar[2]
    mar[2] <- 1
    par(mar=mar)
    par(pty="m")
    plot.new()
    plot.window(xlim=c(0,1), ylim=range(time.vec), xaxs="i",
                yaxs="i")
    rect(0, time.vec[-length(time.vec)], 1, time.vec[-1], col=col.vec)
    axis(4)
    box()
    mar <- mar.orig
    mar[4] <- 1
    par(mar=mar)
    par(pty="s")   
    plot(0:1, 0:1, type="n", main="", xlab=xlab, ylab=ylab)
    rect(0, 0, 1, 1, col=bgcol, lty=0)
  
    for (year in 1:N) {
      time.last <- proc.time()[3]
      time.next <- proc.time()[3]
      while ( (time.next - time.last) < delay ){
        time.next <- proc.time()[3]        
      }
       if (is.finite(intercept[year]) & is.finite(slope[year]))
         abline(intercept[year], slope[year], col=col.vec[year])
       else abline(v=c0[year]/(c0[year]+c1[year]), col=col.vec[year])  
    }

    rect(-0.05, -0.05, 1.05, 0, col="white", lty=0)
    rect(-0.05, -0.05, 0, 1.05, col="white", lty=0)
    rect(-0.05, 1, 1.05, 1.05, col="white", lty=0)
    rect(1, -0.05, 1.05, 1.05, col="white", lty=0)
    box()  
    par(pty="m")
    return(0)
  }
