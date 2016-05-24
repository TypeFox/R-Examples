# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2    
"imagePlotInfo" <- function(..., breaks = NULL, nlevel) {
#NOTE:
# image.plot.info 
# has been renamed as imagePlotInfo to avoid confusion with
# an S3 method
    temp <- list(...)
    #
    xlim <- NA
    ylim <- NA
    zlim <- NA
    poly.grid <- FALSE
    #
    # go through various cases of what these can be
    #
    ##### x,y,z list is first argument
    if (is.list(temp[[1]])) {
        xlim <- range(temp[[1]]$x, na.rm = TRUE)
        ylim <- range(temp[[1]]$y, na.rm = TRUE)
        zlim <- range(temp[[1]]$z, na.rm = TRUE)
        if (is.matrix(temp[[1]]$x) & is.matrix(temp[[1]]$y) & 
            is.matrix(temp[[1]]$z)) {
            poly.grid <- TRUE
        }
    }
    ##### check for polygrid first three arguments should be matrices
    #####
    if (length(temp) >= 3) {
        if (is.matrix(temp[[1]]) & is.matrix(temp[[2]]) & is.matrix(temp[[3]])) {
            poly.grid <- TRUE
        }
    }
    #####  z is passed without an  x and y  (and not a poly.grid!)
    #####
    if (is.matrix(temp[[1]]) & !poly.grid) {
        xlim <- c(0, 1)
        ylim <- c(0, 1)
        zlim <- range(temp[[1]], na.rm = TRUE)
    }
    ##### if x,y,z have all been passed find their ranges.
    ##### holds if poly.grid or not
    #####
    if (length(temp) >= 3) {
        if (is.matrix(temp[[3]])) {
            xlim <- range(temp[[1]], na.rm = TRUE)
            ylim <- range(temp[[2]], na.rm = TRUE)
            zlim <- range(temp[[3]], na.rm = TRUE)           
        }
    }
    # if constant z values perturb the range (1e-8) by epsilon to 
    # avoid other problems in drawing legend later on
    if( !is.na( zlim[1] ) ){
      if( zlim[1] == zlim[2]){
    	if( zlim[1]==0){
    		 zlim[1]<- -1e-8
    		 zlim[2]<- 1e-8}
        else{		 
         delta<- .01*abs(zlim[1])
        zlim[1]<- zlim[1] - delta
        zlim[2]<- zlim[2] + delta
        }
      }
    }
    #### parse x,y,z if they are  named arguments
    # determine if  this is polygon grid (x and y are matrices)
    if (is.matrix(temp$x) & is.matrix(temp$y) & is.matrix(temp$z)) {
        poly.grid <- TRUE
    }
# set limits from the usual $x $y $z format of image object    
    xthere <- match("x", names(temp))
    ythere <- match("y", names(temp))
    zthere <- match("z", names(temp))
    if (!is.na(zthere)) 
        zlim <- range(temp$z, na.rm = TRUE)
    if (!is.na(xthere)) 
        xlim <- range(temp$x, na.rm = TRUE)
    if (!is.na(ythere)) 
        ylim <- range(temp$y, na.rm = TRUE)
        
# overwrite limits with passed values
    if (!is.null(temp$zlim)) 
        zlim <- temp$zlim
    if (!is.null(temp$xlim)) 
        xlim <- temp$xlim
    if (!is.null(temp$ylim)) 
        ylim <- temp$ylim
# At this point xlim, ylim and zlim should be correct 
# using all the different possibilities and defaults for these values
#        
#  Now set up the breaks
    if( is.null(breaks)){
    	midpoints<- seq( zlim[1], zlim[2],,nlevel)
    	delta<- (midpoints[2]- midpoints[1])/2
    	# nlevel +1 breaks with the min and max as midpoints 
    	# of the first and last bins.
    
    	breaks <- c( midpoints[1]- delta, midpoints + delta)
    }        
    list(xlim = xlim, ylim = ylim, zlim = zlim, poly.grid = poly.grid,
       breaks=breaks)
}
# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
# NOTE:
# image.plot.plt<- function(...){
# this function has been renamed as imageplot.setup to avoid confusion with
# an S3 method
#   imageplot.setup(...)}

"imageplot.setup" <- function(x, add = FALSE, legend.shrink = 0.9, 
    legend.width = 1, horizontal = FALSE, legend.mar = NULL, 
    bigplot = NULL, smallplot = NULL, ...) {
    old.par <- par(no.readonly = TRUE)
    if (is.null(smallplot)) 
        stick <- TRUE
    else stick <- FALSE
    if (is.null(legend.mar)) {
        legend.mar <- ifelse(horizontal, 3.1, 5.1)
    }
    # compute how big a text character is
    char.size <- ifelse(horizontal, par()$cin[2]/par()$din[2], 
        par()$cin[1]/par()$din[1])
    # This is how much space to work with based on setting the margins in the
    # high level par command to leave between strip and big plot
    offset <- char.size * ifelse(horizontal, par()$mar[1], par()$mar[4])
    # this is the width of the legned strip itself.
    legend.width <- char.size * legend.width
    # this is room for legend axis labels
    legend.mar <- legend.mar * char.size
    # smallplot is the plotting region for the legend.
    if (is.null(smallplot)) {
        smallplot <- old.par$plt
        if (horizontal) {
            smallplot[3] <- legend.mar
            smallplot[4] <- legend.width + smallplot[3]
            pr <- (smallplot[2] - smallplot[1]) * ((1 - legend.shrink)/2)
            smallplot[1] <- smallplot[1] + pr
            smallplot[2] <- smallplot[2] - pr
        }
        else {
            smallplot[2] <- 1 - legend.mar
            smallplot[1] <- smallplot[2] - legend.width
            pr <- (smallplot[4] - smallplot[3]) * ((1 - legend.shrink)/2)
            smallplot[4] <- smallplot[4] - pr
            smallplot[3] <- smallplot[3] + pr
        }
    }
    if (is.null(bigplot)) {
        bigplot <- old.par$plt
        if (!horizontal) {
            bigplot[2] <- min(bigplot[2], smallplot[1] - offset)
        }
        else {
            bottom.space <- old.par$mar[1] * char.size
            bigplot[3] <- smallplot[4] + offset
        }
    }
    if (stick & (!horizontal)) {
        dp <- smallplot[2] - smallplot[1]
        smallplot[1] <- min(bigplot[2] + offset, smallplot[1])
        smallplot[2] <- smallplot[1] + dp
    }
    return(list(smallplot = smallplot, bigplot = bigplot))
}
# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"crop.image" <- function(obj, loc = NULL, ...) {
    if (is.null(loc)) {
        image.plot(obj, ...)
        loc <- get.rectangle()
    }
    # coerce to midpoints
    m <- nrow(obj$z)
    n <- ncol(obj$z)
    nx <- length(obj$x)
    ny <- length(obj$y)
    if (nx != m) {
        obj$x <- (obj$x[1:m] + obj$x[2:(m + 1)])/2
    }
    if (ny != n) {
        obj$y <- (obj$y[1:n] + obj$x[2:(n + 1)])/2
    }
    # coerce loc to x,y list format if matrix  or data frame
    if (is.matrix(loc) | is.data.frame(loc)) {
        if (ncol(loc) != 2) {
            stop("loc must have two columns\n(for x and y coordinates )")
        }
        loc <- list(x = loc[, 1], y = loc[, 2])
    }
    x <- obj$x
    y <- obj$y
    N <- length(x)
    xr <- range(loc$x)
    xtest <- range(x)
    if (xr[1] < xtest[1] | xr[2] > xtest[2]) {
        stop("cropping outside ranges of x values")
    }
    x1 <- max((1:N)[xr[1] >= x])
    x2 <- min((1:N)[xr[2] <= x])
    N <- length(y)
    yr <- range(loc$y)
    ytest <- range(y)
    if (yr[1] < ytest[1] | yr[2] > ytest[2]) {
        stop("cropping outside ranges of y values")
    }
    y1 <- max((1:N)[yr[1] >= y])
    y2 <- min((1:N)[yr[2] <= y])
    list(x = obj$x[x1:x2], y = obj$y[y1:y2], z = obj$z[x1:x2, 
        y1:y2])
}
average.image <- function(obj, Q = 2) {
    # fast method to sum over a QXQ block in image.
    # Q is the number of elements to average over in each dimension
    # e.g.  Q=5 --  blocks of 25 values are averaged to one grid cell.
    if (is.matrix(obj)) {
        obj <- list(x = 1:nrow(obj), y = 1:ncol(obj), z = obj)
    }
    M <- length(obj$x)
    N <- length(obj$y)
    Mi <- trunc(M/Q)
    Ni <- trunc(N/Q)
    # space to hold results
    z <- matrix(NA, nrow = Mi, ncol = N)
    x2 <- rep(NA, Mi)
    y2 <- rep(NA, Ni)
    indQ <- 1:Q
    # sum over block of rows and handle x grid values
    for (j in 1:Mi) {
        x2[j] <- mean(obj$x[indQ + (j - 1) * Q])
        z[j, ] <- colMeans(obj$z[indQ + (j - 1) * Q, ], na.rm = TRUE)
    }
    # sum over blocks of columns  and average y grid values
    for (k in 1:Ni) {
        y2[k] <- mean(obj$y[indQ + (k - 1) * Q])
        z[, k] <- rowMeans(z[, indQ + (k - 1) * Q], na.rm = TRUE)
    }
    return(list(x = x2, y = y2, z = z[1:Mi, 1:Ni], Q = Q))
}
"get.rectangle" <- function() {
    temp <- locator(2, type = "p", pch = "+")
    rect(temp$x[1], temp$y[1], temp$x[2], temp$y[2])
    temp
}
"half.image" <- function(obj) {
    # coerce to list if a matrix
    if (is.matrix(obj)) {
        obj <- list(x = 1:nrow(obj), y = 1:ncol(obj), z = obj)
    }
    M <- length(obj$x)
    N <- length(obj$y)
    M2 <- trunc(M/2)
    N2 <- trunc(N/2)
    z <- matrix(NA, nrow = M2, ncol = N2)
    ix <- (1:M2) * 2
    iy <- (1:N2) * 2
    x2 <- (obj$x[ix - 1] + obj$x[ix])/2
    y2 <- (obj$y[iy - 1] + obj$y[iy])/2
    return(list(x = x2, y = y2, z = (obj$z[ix - 1, iy] + obj$z[ix - 
        1, iy - 1] + obj$z[ix, iy - 1] + obj$z[ix, iy])/4))
}

pushpin <- function(x, y, z, p.out, height = 0.05, 
    col = "black", text = NULL, adj = -0.1, cex = 1, ...) {
    # project your x,y,z on to the uv plane of the plot
    Sxy1 <- trans3d(x, y, z, p.out)
    Sxy2 <- Sxy1
    hold <- par()$usr
    Sxy2$y <- (hold[4] - hold[3]) * height + Sxy2$y
    # draw the pin
    segments(Sxy1$x, Sxy1$y, Sxy2$x, Sxy2$y, col = "black")
    points(Sxy2, col = col, pch = 19, cex = cex)
    # add a label
    if (!is.null(text)) {
        text(Sxy2$x, Sxy2$y, label = text, adj = adj, cex = cex, 
            ...)
    }
}

designer.colors <- function(n = 256, col = c("darkgreen", 
    "white", "darkred"), x = seq(0, 1,, length(col) ), alpha = 1) {
# generate colors at equal spacings but interpolate to colors at x
    xRange<- range(x)
    xg <- seq(xRange[1], xRange[2],, n)
# convert colors from names  e.g. "magenta" to rgb in [0.1]    
    y.rgb <- t(col2rgb(col))/255
# matrix to hold RGB color values
    temp <- matrix(NA, ncol = 3, nrow = n)
    nColors<- length( col)
    if( nColors != length( x)){
      stop("number of colors needs to be the same as length of x")}
# linear or spline interpolation of RGB color values at x onto xg
    for (k in 1:3) {
        if( nColors > 2){
           hold <- splint(x, y.rgb[, k], xg)}
        else{
           a<-(xRange[2]-xg)/(xRange[2] - xRange[1])
           hold<-  a*y.rgb[1, k] + (1-a)*y.rgb[2, k] }
        # fix up to be in [0,1]
        hold[hold < 0] <- 0
        hold[hold > 1] <- 1
        temp[, k] <- hold
    }
    # convert back to hex
   if(alpha==1){
      return( rgb(temp[, 1], temp[, 2], temp[, 3]))
    }
    else{
      return( rgb(temp[, 1], temp[, 2], temp[, 3], alpha = alpha))
    }
}

#boulder.colors<- c('darkred', 'darkorange',
#                   'white', 'darkgreen', 'darkblue')
"two.colors" <- function(n = 256, start = "darkgreen", 
    end = "red", middle = "white", alpha = 1) {
    designer.colors(n, c(start, middle, end), alpha = alpha)
}

fieldsPlotColors<- function( col, ...){
               N<- length(col)
               image.plot( 1:N, 1, matrix(1:N,N,1), col=col,axes=FALSE, xlab='', ylab='',...)}


imageplot.info<- function (...) 
{
    temp <- list(...)
    xlim <- NA
    ylim <- NA
    zlim <- NA
    poly.grid <- FALSE
    if (is.list(temp[[1]])) {
        xlim <- range(temp[[1]]$x, na.rm = TRUE)
        ylim <- range(temp[[1]]$y, na.rm = TRUE)
        zlim <- range(temp[[1]]$z, na.rm = TRUE)
        if (is.matrix(temp[[1]]$x) & is.matrix(temp[[1]]$y) & 
            is.matrix(temp[[1]]$z)) {
            poly.grid <- TRUE
        }
    }
    if (length(temp) >= 3) {
        if (is.matrix(temp[[1]]) & is.matrix(temp[[2]]) & is.matrix(temp[[3]])) {
            poly.grid <- TRUE
        }
    }
    if (is.matrix(temp[[1]]) & !poly.grid) {
        xlim <- c(0, 1)
        ylim <- c(0, 1)
        zlim <- range(temp[[1]], na.rm = TRUE)
    }
    if (length(temp) >= 3) {
        if (is.matrix(temp[[3]])) {
            xlim <- range(temp[[1]], na.rm = TRUE)
            ylim <- range(temp[[2]], na.rm = TRUE)
            zlim <- range(temp[[3]], na.rm = TRUE)
        }
    }
    if (is.matrix(temp$x) & is.matrix(temp$y) & is.matrix(temp$z)) {
        poly.grid <- TRUE
    }
    xthere <- match("x", names(temp))
    ythere <- match("y", names(temp))
    zthere <- match("z", names(temp))
    if (!is.na(zthere)) 
        zlim <- range(temp$z, na.rm = TRUE)
    if (!is.na(xthere)) 
        xlim <- range(temp$x, na.rm = TRUE)
    if (!is.na(ythere)) 
        ylim <- range(temp$y, na.rm = TRUE)
    if (!is.null(temp$zlim)) 
        zlim <- temp$zlim
    if (!is.null(temp$xlim)) 
        xlim <- temp$xlim
    if (!is.null(temp$ylim)) 
        ylim <- temp$ylim
    list(xlim = xlim, ylim = ylim, zlim = zlim, poly.grid = poly.grid)
}
