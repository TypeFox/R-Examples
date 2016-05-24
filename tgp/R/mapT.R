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


## mapT:
##
## plot the Maximum a Posteriori tree in a tgp-class object,
## or add it to an existing plot -- The proj argument allows
## only some dimensions to be plotted

"mapT" <-
function(out, proj=NULL, slice=NULL, add=FALSE, lwd=2, ...)
{
  ## simple for 1-d data, projection plot
  if(out$d == 1) { proj <- 1; slice <- NULL }
  
  ## otherwise, many options for >= 2-d data

  if(out$d > 2 && !is.null(slice)) { # slice plot
      
    ## will call stop() if something is wrong with the slice
    d <- check.slice(slice, out$d, getlocs(out$X))
    
    ## plot the parts
    tgp.plot.parts.2d(out$parts, d, slice);
    
  } else { # projection plot

    ## will call stop() if something is wrong with the proj
    proj <- check.proj(proj)
    
    ## 1-d projection
    if(length(proj) == 1) {
      if(add == FALSE) plot(out$X[,proj], out$Z, ...)
      tgp.plot.parts.1d(out$parts[,proj], lwd=lwd)

    } else {
    
      ## 2-d projection
      if(add == FALSE) plot(out$X[,proj], ...)
      tgp.plot.parts.2d(out$parts[,proj], lwd=lwd)
    }
  }
}


## tgp.plot.parts.1d:
##
## plot the partitings of 1-d tgp$parts output -- used
## by mapT and plot.tgp

"tgp.plot.parts.1d" <-
function(parts, lwd=2)
{
  j <- 3
  if(is.null(dim(parts))) dp <- length(parts)
  else {
    dp <- nrow(parts)
    parts <- parts[,1]
  }
  is <- seq(2, dp, by=4)
  m <- max(parts[is])
  for(i in is) {
    if(parts[i] == m) next;
    abline(v=parts[i], col=j, lty=j, lwd=lwd);
    j <- j + 1
  }
}


## tgp.plot.parts.2d:
##
## plot the partitings of 2-d tgp$parts output -- used
## by mapT and plot.tgp via tgp.plot.slide and tgp.plot.proj
## the what argument specifies the slice, and trans can make
## rotations

"tgp.plot.parts.2d" <-
function(parts, dx=c(1,2), what=NULL, trans=matrix(c(1,0,0,1), nrow=2),
         col=NULL, lwd=3)
{
  if(length(what) > 0) {
    indices <- c()
    for(i in seq(1,nrow(parts),4)) {
      opl <- i+2; opr <- i+3;
      if(parts[opl,what$x] == 104 && parts[opr,what$x] == 102
         && what$z >= parts[i,what$x] && what$z <= parts[i+1,what$x]) {
        indices <- c(i, indices)
      } else if(parts[opl,what$x] == 105 && parts[opr,what$x] == 102
                && what$z > parts[i,what$x] && what$z <= parts[i+1,what$x]) {
        indices <- c(i, indices)
      } 
    }
  } else {
    indices <- seq(1,dim(parts)[1],4);
  }
  
  j <- 1
  for(i in indices) {
    a <- parts[i,dx[1]]; b <- parts[i+1,dx[1]];
    c <- parts[i,dx[2]]; d <- parts[i+1,dx[2]];
    x <- c(a, b, b, a, a);
    y <- c(c, c, d, d, c);
    xy <- as.matrix(cbind(x,y)) %*% trans
    if(is.null(col)) { lines(xy, col=j, lty=j, lwd=lwd); }
    else { lines(xy, col=col, lty=1, lwd=lwd); }
    j <- j+1
  }
}
