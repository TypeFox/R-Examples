
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU Library General Public License for more details.
#
# You should have received A copy of the GNU Library General 
# Public License along with this library; if not, write to the 
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
# MA  02111-1307  USA

# Copyrights (C)
# for this R-port: 
#   1999 - 2007, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:                 DESCRIPTION:
#  squareBinning             Square binning of irregularly spaced points
#  plot                      S3 Method for plotting square binned points
# FUNCTION:                 DESCRIPTION:
#  hexBinning                Hexagonal binning of irregularly spaced points
#  plot                      S3 Method for plotting hexagonal binned points
################################################################################


################################################################################
# FUNCTION:                 DESCRIPTION:
#  squareBinning             Square binning of irregularly spaced points
#  plot                      S3 Method for plotting square binned points


squareBinning = 
function(x, y = NULL, bins = 30)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns 2D Histogram Counts
    
    # Arguments:
    #   x, y - two vectors of coordinates of data. If y is NULL then x 
    #       is assumed to be a two column matrix, where the first column 
    #       contains the x data, and the second column the y data. 
    #       'timeSeries' objects are also allowed as input.
    #   bins - number of bins in each dimension, may be a scalar or a 2
    #       element vector. The default value is 20.
    
    # Value:
    #    A list with three elements x, y, and z. x and y are vectors
    #       spanning the two dimensioanl grid and z the corresponding
    #       matrix. The output can directly serve as input to the
    #       plotting functions image, contour and persp.
    
    # Example:
    #   sB = squareBinning(x = rnorm(1000), y = rnorm(1000)); plot(sB)
   
    # Note:
    #   Partly copied from R Package gregmisc, function 'hist2d'.
    
    # FUNCTION:
    
    # 2D Histogram Counts:
    if (is.null(y)) {
        x = as.matrix(x)
        y = x[, 2]
        x = x[, 1]
    } else {
        x = as.vector(x)
        y = as.vector(y)
    }
    data = cbind(x, y)
    
    # Bins:
    n = bins
    if (length(n) == 1) {
        nbins = c(n, n)
    } else {
        nbins = n
    }
    
    # Binning:
    xo = seq(min(x), max(x), length = nbins[1])
    yo = seq(min(y), max(y), length = nbins[2])
    xvals = xo[-1] - diff(xo)/2
    yvals = yo[-1] - diff(yo)/2
    ix = findInterval(x, xo)
    iy = findInterval(y, yo)
    xcm = ycm = zvals = matrix(0, nrow = nbins[1], ncol = nbins[2])
    
    for (i in 1:length(x)) {
        zvals[ix[i], iy[i]] = zvals[ix[i], iy[i]] + 1
        xcm[ix[i], iy[i]] = xcm[ix[i], iy[i]] + x[i]
        ycm[ix[i], iy[i]] = ycm[ix[i], iy[i]] + y[i]
    }

    # Reduce to non-empty cells:
    u = v = w = ucm = vcm = rep(0, times = nbins[1]*nbins[2])
    L = 0
    for (i in 1:(nbins[1]-1)) {
        for (j in 1:(nbins[2]-1)) {
            if (zvals[i, j] > 0) {
                L = L + 1
                u[L] = xvals[i]
                v[L] = yvals[j]
                w[L] = zvals[i, j] 
                ucm[L] = xcm[i, j]/w[L]
                vcm[L] = ycm[i, j]/w[L]
            }    
        }    
    }
    length(u) = length(v) = length(w) = L    
    length(ucm) = length(vcm) = L    
    
    ans = list(x = u, y = v, z = w, xcm = ucm, ycm = vcm, bins = bins,
        data = data)
    class(ans) = "squareBinning"
    
    # Return Value:
    ans
    
}


# ------------------------------------------------------------------------------


plot.squareBinning =
function(x, col = heat.colors(12), addPoints = TRUE, addRug = TRUE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Plot square binned data points
    
    # FUNCTION:
    
    # Binning:
    X = x$x
    Y = x$y
    
    # Plot Center Points:
    plot(X, Y, type = "n", ...)
      
    # Create Hexagon Coordinates:
    rx = min(diff(unique(sort(X))))/2
    ry = min(diff(unique(sort(Y))))/2
    u = c(-rx, rx,  rx, -rx)
    v = c( ry, ry, -ry, -ry)

    # Create Color Palette:
    N = length(col)
    Z = x$z
    zMin = min(Z)
    zMax = max(Z)
    Z = (Z - zMin)/(zMax - zMin)
    Z = trunc(Z*(N-1)+1)
    
    # Add Colored Hexagon Polygons:
    for (i in 1:length(X)) {
        polygon(u+X[i], v+Y[i], col = col[Z[i]], border = "white")
    }
    
    # Add Center of Mass Points:
    if (addPoints) {
        points(x$xcm, x$ycm, pch = 19, cex = 1/3, col = "black")
    }
    
    # Add rug:
    if (addRug) {
        rug(x$data[, 1], ticksize = 0.01, side = 3)
        rug(x$data[, 2], ticksize = 0.01, side = 4)
    }
    
    
    # Return Value:
    invisible(NULL)
}


################################################################################
# FUNCTION:                 DESCRIPTION:
#  hexBinning                Hexagonal binning of irregularly spaced points
#  plot                      S3 Method for plotting hexagonal binned points


hexBinning = 
function(x, y = NULL, bins = 30)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Does a hexagonal binning of data points
    
    # Arguments:
    #   x, y - two vectors of coordinates of data. If y is NULL then x 
    #       is assumed to be a two column matrix, where the first column 
    #       contains the x data, and the second column the y data. 
    #       'timeSeries' objects are also allowed as input.
    #   bins - number of bins in each dimension, may be a scalar or a 2
    #       element vector. The default value is 20.
    
    # Example:
    #   hB = hexBinning(x = rnorm(10000), y = rnorm(10000)); plot(hB)
       
    # FUNCTION:
    
    # Extract Series:
    if (is.null(y)) {
        x = as.matrix(x)
        y = x[, 2]
        x = x[, 1]
    } else {
        x = as.vector(x)
        y = as.vector(y)
    }
    data = cbind(x, y)
    
    # Set Parameters:
    shape = 1
    n = length(x)
    xbnds = range(x)
    ybnds = range(y)
    jmax = floor(bins + 1.5001)
    c1 = 2 * floor((bins *shape)/sqrt(3) + 1.5001)
    imax = trunc((jmax*c1 -1)/jmax + 1)
    lmax = jmax * imax     
    cell = cnt = xcm = ycm = rep(0, times = max(n, lmax))
    xmin = xbnds[1]
    ymin = ybnds[1]
    xr = xbnds[2] - xmin
    yr = ybnds[2] - ymin
    c1 = bins/xr
    c2 = bins*shape/(yr*sqrt(3.0))
    jinc = jmax
    lat = jinc + 1
    iinc = 2*jinc
    con1 = 0.25
    con2 = 1.0/3.0

    # Count Bins:
    for ( i in 1:n ) {
        sx = c1 * (x[i] - xmin)
        sy = c2 * (y[i] - ymin)
        j1 = floor(sx + 0.5)
        i1 = floor(sy + 0.5)
        dist1 = (sx-j1)^2 + 3.0*(sy-i1)^2
        if( dist1 < con1) {
            L = i1*iinc + j1 + 1
        } else if (dist1 > con2) {
            L = floor(sy)*iinc + floor(sx) + lat
        } else {
            j2 = floor(sx)
            i2 = floor(sy)
            test = (sx-j2 -0.5)^2 + 3.0*(sy-i2-0.5)^2
            if ( dist1 <= test ) {
                L = i1*iinc + j1 + 1
            } else {
                L = i2*iinc + j2 + lat
            }
        }
        cnt[L] = cnt[L]+1
        xcm[L] = xcm[L] + (x[i] - xcm[L])/cnt[L]
        ycm[L] = ycm[L] + (y[i] - ycm[L])/cnt[L]
    }

    # Reduce to Non-Empty Cells:
    nc = 0
    for ( L in 1:lmax ) {
        if(cnt[L] > 0) {
            nc = nc + 1
            cell[nc] = L
            cnt[nc] = cnt[L]
            xcm[nc] = xcm[L]
            ycm[nc] = ycm[L]
        }
    }
    bnd = c(imax, jmax)
    bnd[1] = (cell[nc]-1)/bnd[2] + 1
    length(cell) = nc
    length(cnt) = nc
    length(xcm) = nc
    length(ycm) = nc
    if(sum(cnt) != n) warning("Lost counts in binning")
    
    # Compute Positions:
    c3 = diff(xbnds)/bins
    ybnds = ybnds
    c4 = (diff(ybnds) * sqrt(3))/(2 * shape * bins)
    cell = cell - 1
    i = cell %/% jmax
    j = cell %% jmax
    y = c4 * i + ybnds[1]
    x = c3 * ifelse(i %% 2 == 0, j, j + 0.5) + xbnds[1]
    
    # Result:
    ans = list(x = x, y = y, z = cnt, xcm = xcm, ycm = ycm, bins = bins,
        data = data)
    class(ans) = "hexBinning"
    
    # Return Value:
    ans
   
}


# ------------------------------------------------------------------------------


plot.hexBinning =
function(x, col = heat.colors(12), addPoints = TRUE, addRug = TRUE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Plot hexagonal binned data points
    
    # Example:
    #   hexPlot(rnorm(1000), rnorm(1000), bins = 20)
    
    # FUNCTION:
    
    # Binning:
    X = x$x
    Y = x$y
    
    # Plot Center Points:
    plot(X, Y, type = "n", ...)
      
    # Create Hexagon Coordinates:
    rx = min(diff(unique(sort(X))))
    ry = min(diff(unique(sort(Y))))
    rt = 2*ry
    u = c(rx,  0, -rx, -rx,   0,  rx)
    v = c(ry, rt,  ry, -ry, -rt, -ry) / 3

    # Create Color Palette:
    N = length(col)
    z = x$z
    zMin = min(z)
    zMax = max(z)
    Z = (z - zMin)/(zMax - zMin)
    Z = trunc(Z*(N-1)+1)
    
    # Add Colored Hexagon Polygons:
    for (i in 1:length(X)) {
        polygon(u+X[i], v+Y[i], col = col[Z[i]], border = "white")
    }
    
    # Add Center of Mass Points:
    if (addPoints) {
        points(x$xcm, x$ycm, pch = 19, cex = 1/3, col = "black")
    }
    
    # Add rug:
    if (addRug) {
        rug(x$data[, 1], ticksize = 0.01, side = 3)
        rug(x$data[, 2], ticksize = 0.01, side = 4)
    }
    
    # Return Value:
    invisible(NULL)
}


################################################################################


