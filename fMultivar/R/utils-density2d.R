
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:             DESCRIPTION:
#  density2d             Returns 2D Kernel Density Estimates
#  hist2d                Returns 2D Histogram Counts
################################################################################


density2d <- 
    function (x, y = NULL, n = 20, h = NULL, limits = c(range(x), range(y))) 
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns 2D Kernel Density Estimates
    
    # Arguments:
    #   x, y - two vectors of coordinates of data. If y is NULL then x
    #       is assumed to be a two column matrix, where the first column 
    #       contains the x data, and the second column the y data. 
    #   n - Number of grid points in each direction. 
    #   h - a vector of bandwidths for x and y directions. Defaults to
    #       normal reference bandwidth. 
    #   limits - the limits of the rectangle covered by the grid.    
    
    # Value:
    #    A list with three elements x, y, and z. x and y are vectors
    #       spanning the two dimensioanl grid and z the corresponding
    #       matrix. The output can directly serve as input to the
    #       plotting functions image, contour and persp.
    
    # Details:
    #   Two-dimensional kernel density estimation with an axis-aligned
    #   bivariate normal kernel, evaluated on a square grid.
    
    # Note:
    #   Partly copied from R Package MASS, function 'kde2d'.
    
    # Reference:
    #   Venables, W.N., Ripley, B. D. (2002); 
    #       Modern Applied Statistics with S.
    #       Fourth edition, Springer.
    
    # FUNCTION:
    
    # Settings:
    lims = limits
    if (is.null(y)) {
        y = x[, 2]
        x = x[, 1]
    }
    
    # Bandwidth:
    .bandwidth.nrd = function (x) {
        r = quantile(x, c(0.25, 0.75))
        h = (r[2] - r[1])/1.34
        4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5) }
        
    # Kernel Density Estimator:
    nx = length(x)
    if (length(y) != nx) stop("Data vectors must be the same length")
    gx = seq(lims[1], lims[2], length = n)
    gy = seq(lims[3], lims[4], length = n)
    if (is.null(h)) h = c(.bandwidth.nrd(x), .bandwidth.nrd(y))
    h = h/4
    ax = outer(gx, x, "-")/h[1]
    ay = outer(gy, y, "-")/h[2]
    z = matrix(dnorm(ax), n, nx) %*% t(matrix(dnorm(ay), n, 
        nx))/(nx * h[1] * h[2])
    
    # Return Value:    
    list(x = gx, y = gy, z = z)
}


# ------------------------------------------------------------------------------


hist2d <-  
    function(x, y = NULL, n = c(20, 20))
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns 2D Histogram Counts
    
    # Arguments:
    #   x, y - two vectors of coordinates of data. If y is NULL then x 
    #       is assumed to be a two column matrix, where the first column 
    #       contains the x data, and the second column the y data. 
    #   n - number of bins in each dimension, may be a scalar or a 2
    #       element vector. The default value is 20.
    
    # Value:
    #    A list with three elements x, y, and z. x and y are vectors
    #       spanning the two dimensioanl grid and z the corresponding
    #       matrix. The output can directly serve as input to the
    #       plotting functions image, contour and persp.
   
    # Note:
    #   Partly copied from R Package gregmisc, function 'hist2d'.
    
    # FUNCTION:
    
    # 2D Histogram Counts:
    if (is.null(y)) {
        y = x[, 2]
        x = x[, 1]
    }
    if (length(n) == 1) {
        nbins = c(n, n)
    } else {
        nbins = n
    }
    nas = is.na(x) | is.na(y)
    x.cuts = seq(from = min(x, y), to = max(x,y), length = nbins[1]+1)
    y.cuts = seq(from = min(x, y), to = max(x,y), length = nbins[2]+1)
    index.x = cut(x, x.cuts, include.lowest = TRUE)
    index.y = cut(y, y.cuts, include.lowest = TRUE)
    m = matrix(0, nrow=nbins[1], ncol = nbins[2],
        dimnames = list( levels(index.x), levels(index.y) ) )
    for ( i in 1:length(index.x) ) {
        m[index.x[i], index.y[i] ] = m[index.x[i], index.y[i] ] + 1
    }
    xvals = x.cuts[1:nbins[1]]
    yvals = y.cuts[1:nbins[2]]

    # Return Value:
    list(x = xvals, y = yvals, z = m)
}


################################################################################

