#### Author: Adriano Zanin Zambom
#### contact: adriano.zambom@gmail.com
#### last modified: 21/Set/2015
####
#### papers of reference: 
#### Zambom, A. Z. and Akritas, M. G. (2012). a) Nonparametric Model Checking and Variable Selection. Statistica Sinica, v. 24, pp. 1837
#### Zambom, A. Z. and Akritas, M. G. (2012). b) Signicance Testing and Group Variable Selection. Journal of Multivariate Analysis, v. 133, pp. 51

 

plot3d.localpoly.reg <- function(X,Y, bandwidth = "CV", gridsize = 30, degree.pol = 0, kernel.type = "epanech", gridsurface = 30, xlab=expression(X_1), ylab=expression(X_2), zlab=expression(Y), theta = 30, phi = 30, expand = 0.5, col = "lightblue",
ltheta = 120, shade = 0.75, ticktype = "detailed", pch = 16, ... )
{
    if (is.null(dim(X)))
        stop("\n\nplot3d for localpoly.reg is only available for X with exactly 2 dimensions") else
    if (dim(X)[2] != 2)
        stop("\n\nplot3d for localpoly.reg is only available for X with exactly 2 dimensions") else
    if (!is.null(dim(bandwidth)))
       if(length(bandwidth) != gridsurface^2*2)
          stop("\n\n if inputing a matrix of bandwidths, it must be of dimensions gridsurface^2 by 2")
    
    
        n = dim(X)[1]
        rangeX = max(X[,1]) - min(X[,1])
        rangeY = max(X[,2]) - min(X[,2])


        
        W = matrix(0,gridsurface^2,2)         ## the rows of W are the points equally spaced at a grid which 
        for (i in 1:gridsurface)              ## starts at minimum of dist(X) and goes up until max dist(X)
           for (j in 1:gridsurface)
              W[(i-1)*gridsurface + j,] = c(min(X[,1])+i*rangeX/gridsurface,min(X[,2])+j*rangeY/gridsurface)
        
        
        m1 = localpoly.reg(X,Y, points=W, bandwidth=bandwidth, gridsize = gridsize, degree.pol=degree.pol, kernel.type=kernel.type)

        
        z <- matrix(NA, gridsurface, gridsurface)   ## z is a matrix that needs to be passed to persp, containing
                                              ## the estimated points of each grid point (grid from W)
        for (i in 1:gridsurface)
           for (j in 1:gridsurface)
              if (m1$predicted[(i-1)*gridsurface + j] != 0)
                 z[i,j] = m1$predicted[(i-1)*gridsurface + j];
        
        
        
        op <- par(bg = "white")
        persp(W[seq(1,length(W[,1]),by=gridsurface),1],W[1:gridsurface,2], z, xlim = range(X[,1]), ylim = range(X[,2]), zlim = range(Y), xlab=xlab, ylab=ylab, zlab=zlab, theta = theta, phi = phi, expand = expand, col = col, 
        ltheta = ltheta, shade = shade, ticktype = ticktype, ...)-> res
        
        for (i in 1:n)
           points(trans3d(X[i,1], X[i,2], Y[i], pmat = res), col = 2, pch = pch)
    
        result = list()  
        result$x = X
        result$y = Y
        result$points = W
        result$bandwidth = m1$bandwidth
        result$predicted = z

        return(result)
    
    
}

