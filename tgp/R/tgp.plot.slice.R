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


## tgp.plot.slice:
##
## plot the 2-d slice of the tgp-class object out specified by the slice
## argument, and other usual plotting arguments and specified center
## and as (error) specifications

"tgp.plot.slice" <-
function(out, pparts=TRUE, slice=NULL, map=NULL, as=NULL, center="mean",
         layout="both", main=NULL, xlab=NULL, ylab=NULL, zlab=NULL,
         pc="pc", gridlen=40, span=0.1,...)
{
  ## choose center as median or mean (i.e., X & Z data)
  ## (this hasn't been tested since the addition of the tgp.choose.center() function
  center <- tgp.choose.center(out, center);
  Z.mean <- center$Z
  cname <- center$name;
  X <- center$X
   
  ## get X locations for calculating slice
  locs <- getlocs(X)

  ## will call stop() if something is wrong with the slice
  d <- check.slice(slice, out$d, locs)
  
  ## deal with axis labels
  if(is.null(xlab)) xlab <- names(out$X)[d[1]]
  if(is.null(ylab)) ylab <- names(out$X)[d[2]]
  if(is.null(zlab)) zlab <- out$response
  fixed <- names(out$X)[slice$x]; to <- slice$z
  slice.str <- paste("(", fixed, ") fixed to (", to, ")", sep="")
  smain <- paste(main, " ", zlab, " ", cname, ", with ", slice.str, sep="")
  
  ## for ALC and EGO plotting
  as <- tgp.choose.as(out, as);
  XX <- as$X
  ZZ.q <- as$criteria
  emain <- paste(main, " ", zlab, " ", as$name, ", with ",  slice.str, sep="")
  ##emain <- paste(main, zlab, as$name)

  ## depict the slice in terms of index variables p*
  if(length(slice$x) > 1) {
    p <- seq(1,nrow(X))[apply(X[,slice$x] == slice$z, 1, prod) == 1]
    pp <- seq(1,nrow(XX))[apply(XX[,slice$x] == slice$z, 1, prod) == 1]
    pn <- seq(1,out$n)[apply(out$X[,slice$x] == slice$z, 1, prod) == 1]
    ppn <- seq(1,out$nn)[apply(out$XX[,slice$x] == slice$z, 1, prod) == 1]
  } else {
    ppn <- seq(1,out$nn)[(out$XX[,slice$x] == slice$z)]
    pn <- seq(1,out$n)[out$X[,slice$x] == slice$z]
    p <- seq(1,nrow(X))[X[,slice$x] == slice$z]
    pp <- seq(1,nrow(XX))[XX[,slice$x] == slice$z]
  }
  
  ## check to makes sure there is actually some data in the slice
  if(length(p) == 0) {
    print(slice)
    stop("no points in the specified slice\n")
  }
  
  ## prepare for plotting
  if(layout == "both") par(mfrow=c(1,2), bty="n")
  ## else par(mfrow=c(1,1), bty="n")
    
  Xd.1 <- X[,d[1]]; Xd.2 <- X[,d[2]]
  XXd.1 <- XX[,d[1]]; XXd.2 <- XX[,d[2]]
  
  if(pc == "c") { # double-image plot
    if(layout == "both" || layout == "surf") {
      slice.image(Xd.1,Xd.2,p,Z.mean,main=smain,xlab=xlab,ylab=ylab,
                  gridlen=gridlen,span=span,...)
      if(pparts & !is.null(out$parts)) { tgp.plot.parts.2d(out$parts, d, slice); }
      if(length(pn) > 0) points(out$X[pn,d[1]], out$X[pn,d[2]], pch=20)
      if(length(ppn) > 0) points(out$XX[ppn,d[1]], out$X[ppn,d[2]], pch=21)
    }
    if(layout == "both" || layout == "as") {
      slice.image(XXd.1,XXd.2,pp,ZZ.q,main=emain,xlab=xlab,ylab=ylab,
                      gridlen=gridlen,span=span,...)
      if(pparts & !is.null(out$parts)) { tgp.plot.parts.2d(out$parts, d, slice); }
      if(length(pn) > 0) points(out$X[pn,d[1]], out$X[pn,d[2]], pch=20)
      if(length(ppn) > 0) points(out$XX[ppn,d[1]], out$XX[ppn,d[2]], pch=21)
      if(substr(as$name,1,1) == "I")
        text(out$XX[ppn,d[1]], out$XX[ppn,d[2]], labels=out$improv[ppn,2], ...)
    }
  } else if(pc == "pc") {	# perspactive and image plot
    if(layout == "both" || layout == "surf")
      slice.persp(Xd.1,Xd.2,p,Z.mean,main=smain,xlab=xlab,ylab=ylab,zlab=zlab,
                  gridlen=gridlen,span=span,...)
    if(layout == "both" || layout == "as") {
      slice.image(XXd.1,XXd.2,pp,ZZ.q,main=emain,xlab=xlab,ylab=ylab,
                  gridlen=gridlen,span=span,...)
      if(length(pn) > 0) points(out$X[pn,d[1]], out$X[pn,d[2]], pch=20)
      if(length(ppn) > 0) points(out$XX[ppn,d[1]], out$XX[ppn,d[2]], pch=21)
      if(pparts & !is.null(out$parts)) { tgp.plot.parts.2d(out$parts, d, slice); }
      if(substr(as$name,1,1) == "I")
        text(out$XX[,proj[1]], out$XX[,proj[2]], labels=out$improv[ppn,2], ...)
    }
  }
}


## slice.contour:
##
## contour plot of the slice or projection -- not currently
## used in any tgp plotting function

"slice.contour" <-
function(x,y,p,z,levels=NULL,xlab="x",ylab="y",main="",xlim=NULL,ylim=NULL, ...)
{
  g <- slice.interp(x,y,p,z,xlim,ylim)
  if(missing(ylim)) ylim <- range(y)
  if(missing(xlim)) xlim <- range(x)
  if(is.null(levels)) {
    contour(g,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,...)
  } else {
    contour(g,levels=levels,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,...)
  }
}


## slice.image:
##
## image plot of the slice or projection -- used in pc="c"

"slice.image" <-
function(x,y,p,z,xlim=NULL, ylim=NULL,
         gridlen=c(40,40), span=0.05, col=terrain.colors(128), ...)
{
  g <- slice.interp(x,y,p,z,xlim,ylim,gridlen=gridlen,span=span)
  if(missing(ylim)) ylim <- range(y)
  if(missing(xlim)) xlim <- range(x)
  image(g, col=col,xlim=xlim,ylim=ylim,...)
}


## slice.image.contour:
##
## double image and contour plot of the slice or projection -- not
## currently used in any tgp plotting function

"slice.image.contour" <-
function(x,y,p,z, xlim=NULL, ylim=NULL,
         gridlen=c(40,40), span=0.05, ...)
{
  g <- slice.interp(x,y,p,z,xlim,ylim,gridlen=gridlen,span=span)
  if(missing(ylim)) ylim <- range(y)
  if(missing(xlim)) xlim <- range(x)
  image(g, col=terrain.colors(128),xlim=xlim,ylim=ylim,...)
  contour(g, add=TRUE,...)
}


## slice.persp:
##
## perspective plot of the slice or projections -- used when
## pc="p"

"slice.persp" <-
function(x,y,p,z,theta=-30,phi=20,xlim=NULL, ylim=NULL,
         gridlen=c(40,40), span=0.05, ...)
{
  g <- slice.interp(x,y,p,z,xlim,ylim,gridlen=gridlen,span=span)
  if(missing(ylim)) ylim <- range(y)
  if(missing(xlim)) xlim <- range(x)
  persp(g, theta=theta, phi=phi, axes=TRUE, box=TRUE, xlim=xlim, ylim=ylim, ...)
}


## slice.interp:
##
## interpolate the x, y, z data specified onto a regular 2-d
## grid, perhaps making a slice specified by the p-vector indicating
## which entries of x, y, and z should be used.  This is necessary
## in order to plot using persp, contour, image, etc.  
## loess is used for interpolation

"slice.interp" <-
function(x, y, p=NULL, z, xlim=NULL, ylim=NULL, gridlen=c(40,40),
         span=0.05, ...)
{
  ## check gridlen
  if(length(gridlen) != 2) stop("length(gridlen) should be 2")
  
  # check and/or default the projection parameter p
  if(is.null(p)) p <- 1:length(x)
  else p <- as.integer(p)
  if(any(p <= 0) || any(p > length(x)))
    stop("invalid p (third arg: value unknown)")

  # make projection
  x <- x[p]; y <- y[p]; z <- z[p]
  if(!is.null(xlim)) { # crop (zoom in) x
    p <- x>=xlim[1] & x<=xlim[2]
    x <- x[p]; y <- y[p]; z <- z[p]
  }
  if(!is.null(ylim)) { # crop (zoom in) y
    p <- y>=ylim[1] & y<=ylim[2]
    x <- x[p]; y <- y[p]; z <- z[p]
  }
  
  # use loess
  return(interp.loess(x,y,z, duplicate="mean", gridlen=gridlen, span=span, ...))
}


## check.slice:
##
## checks to make sure the slice argument to plot.tgp is of a
## format that make sens for the input dimension and data locations
## provided from the getlocs function

"check.slice" <- 
function(slice, dim, locs)
{
  ## check to make sure the slice requested is valid
  numfix <- dim-2;
  if(length(slice$x) != numfix && length(slice$x) == length(slice$z)) {
    print(locs)
    stop(paste("must fix", numfix, "variables, each at one of the above locations\n"))
  }

  ## check to make sure enough dimensions have been fixed
  d <- setdiff(seq(1:dim), slice$x)
  if(length(d) != 2) 
    stop(paste(length(d)-2, "more dimensions need to be fixed\n", sep=""))

  ## will stop if the slice is not ok,
  ## otherwise returns the remaining (unfixed) dimensions
  return(d)
}

## getlocs:
##
## get the grid of locations for the data -- these are the locations
## used in the locs argument of check.slice

"getlocs" <- 
function(X)
{
  db <- dim(X);
  Xsort <- apply(X, 2, sort)
  unique <- (Xsort[1:(db[1]-1),] != Xsort[2:db[1],])
  locs.list <- list()
  for(i in 1:db[2]) {
    locs <- c(Xsort[unique[,i],i], Xsort[db[1],i])
    count <- rep(0,length(locs))
    for(j in 1:length(locs)) {
      count[j] = sum(Xsort[,i] == locs[j])
    }
    ll.i <- list(locs=locs,count=count)
    locs.list[[i]] <- ll.i
  }
  return(locs.list)
}


## interp.loess:
##
## interolate x,y,z onto a regular 2-d grid of size gridlen.
## this function is meant to mimic the interp function in the
## akima library which can be buggy.  It luse a loess smoother
## instead, with the span provided

"interp.loess" <-
function(x, y, z, gridlen = c(40,40), span=0.1, ...)
{
  ## check the gridlen argument
  if(length(gridlen) == 1) gridlen <- rep(gridlen, 2)
  if(length(gridlen) != 2) stop("length(gridlen) should be 2")

  if(length(x) != length(y) && length(y) != length(z))
    stop("length of x, y and z must be equal")

  if(length(x) < 30 && span < 0.5) {
    warning("with less than 30 points, suggest span >> 0.5 or use akima",
            immediate. = TRUE)
    cat(paste("currently trying span =", span, "for", length(x), "points\n"))
  }
  
  xo <- seq(min(x), max(x), length=gridlen[1])
  yo <- seq(min(y), max(y), length=gridlen[2])

  xyz.loess <-
    suppressWarnings(loess(z ~ x + y, data.frame(x=x, y=y), span=span, ...))

  g <- expand.grid(x=xo, y=yo)
  g.pred <- predict(xyz.loess, g)
  return(list(x=xo, y=yo, z=matrix(g.pred, nrow=gridlen)))
}
