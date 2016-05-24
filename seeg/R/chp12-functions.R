Okriging <- function (dataset, vario, step, maxdist,border.sw=FALSE,border.poly="none"){

# dataset: columns 1 and 2 are x,y coordinates, 3 is variable
# vario is variogram model, step is interval for the prediction grid
# maxdist is max distance for prediction

# extract names of coord and variable
 x <- names(dataset)[1]
 y <- names(dataset)[2] 
 v <- names(dataset)[3] 
 
# First, select a grid for the prediction.
# Use min and max of the original dataset and a distance step

grid <- list(x=seq(min(dataset[,1]),max(dataset[,1]),by=step),
             y=seq(min(dataset[,2]),max(dataset[,2]),by=step))

# Get the range and span of x for this prediction grid 
 
grid$xr <- range(grid$x)
grid$xs <- grid$xr[2] - grid$xr[1]

#Get the range and span of y for this prediction grid 

grid$yr <- range(grid$y)
grid$ys <- grid$yr[2] - grid$yr[1]

# Then build two matrices to arrange grid coordinates for prediction 

xmat <- matrix(grid$x, length(grid$x), length(grid$y))
ymat <- matrix(grid$y, length(grid$x), length(grid$y), byrow=TRUE)

#Reconvert to vectors, bind columns and convert to dataframe

grid$xy <- data.frame(cbind(c(xmat), c(ymat)))
colnames(grid$xy) <- c("x", "y")

# convert grid and datset to point patterns

grid$point <- point(grid$xy)
data.point <- point(dataset,x=x,y=y)

# Apply function krige( ) to predict over this grid point pattern using measured points in the dataset, 
# using the variogram model and over distances from measured points not to exceed maxdist.
if(border.sw==FALSE)grid$krige <- krige(grid$point,data.point,v,vario, maxdist=maxdist,extrap=FALSE)
else grid$krige <- krige(grid$point,data.point,v,vario, maxdist=maxdist,extrap=FALSE,border=border.poly)

result <- data.frame(x=grid$krige$x, y=grid$krige$y, zhat=grid$krige$zhat, varhat=grid$krige$sigma2hat)

return(result)
}

plotkriged <- function (dataset, kriged, outpdf="dataset-kriged.pdf",border.sw=FALSE,border.poly="none"){

# dataset: columns 1 and 2 are x,y coordinates
x <- names(dataset)[1]
y <- names(dataset)[2]

# prediction grid given by kriged dataset 
grid <- list(x=unique(kriged$x),y=unique(kriged$y))

# Get the range and span of x and y for the prediction grid 
grid$xr <- range(grid$x)
grid$xs <- grid$xr[2] - grid$xr[1]
grid$yr <- range(grid$y)
grid$ys <- grid$yr[2] - grid$yr[1]

# Build two matrices to arrange grid coordinates for prediction 
xmat <- matrix(grid$x, length(grid$x), length(grid$y))
ymat <- matrix(grid$y, length(grid$x), length(grid$y), byrow=TRUE)

#Reconvert to vectors, bind columns and convert to dataframe
grid$xy <- data.frame(cbind(c(xmat), c(ymat)))
colnames(grid$xy) <- c("x", "y")

# Convert grid and dataset to point patterns
grid$point <- point(grid$xy)
data.point <- point(dataset,x=x,y=y)

pdf(file=outpdf,7,7)
# Declare a square plot, and find maximum to limit the plot.
par(pty="s")
grid$max <- max(grid$xs, grid$ys)
plot(grid$xy, type="n", xlim=c(grid$xr[1], grid$xr[1]+grid$max),
                    ylim=c(grid$yr[1], grid$yr[1]+grid$max))
# Form matrix of kriged values for the image plot
Zest <- matrix(kriged$zhat, length(grid$x),length(grid$y))

# Plot using image and then overlay a contour line map
image(grid$x, grid$y, Zest, col=grey(seq(0.5,1,0.02)),add=TRUE)
contour(grid$x, grid$y, Zest, col=1,add=TRUE)

# Plot of the original point pattern (measured points) 
points(data.point)
title("Kriged values")
if(border.sw==TRUE) lines(border.poly$x,border.poly$y,col=1,lwd=2.5)
# plot variance
par(pty="s")
plot(grid$xy, type="n", xlim=c(grid$xr[1], grid$xr[1]+grid$max),
                    ylim=c(grid$yr[1], grid$yr[1]+grid$max))
Sigmaest <- matrix(kriged$varhat,length(grid$x),length(grid$y))
image(grid$x,grid$y, Sigmaest, col=grey(seq(0.5,1,0.02)),add=TRUE)
contour(grid$x,grid$y, Sigmaest, col=1, add=TRUE)
points(data.point)
title("Variance Kriging error")
if(border.sw==TRUE) lines(border.poly$x,border.poly$y,col=1,lwd=2.5)
dev.off()

}

scan.map.ras <- function(filename){
  # function to scan semi-geoeas files in raster formats
  ncols <- as.numeric(scan(filename, skip=1, what=c("",0), nlines=1)[2])
  nrows <- as.numeric(scan(filename, skip=2, what=c("",0), nlines=1)[2])
  size <- as.numeric(scan(filename, skip=3, what=c("",0), nlines=1)[2])
  nvar <-  as.numeric(scan(filename, skip=4, what=c("",0), nlines=1)[2])
  namesvar <- as.character(structure(1:nvar))
  variab.ras <- list()
  for(i in 1:nvar)
  namesvar[i] <- scan(filename, skip=5+(i-1), what=c("",""), nlines=1)[2]
  variab <- matrix(scan(filename, skip=5+nvar), ncol=nvar, byrow=TRUE)
   x <- seq(0+size/2,size*ncols-size/2,by=size)
   y <- seq(0+size/2,size*nrows-size/2,by=size)

  for(vcol in 1:nvar){
   variab.ras[[vcol]] <- matrix(variab[,vcol], ncol=ncols,byrow=TRUE)
   z <- matrix(t(variab.ras[[vcol]]),length(x),length(y))
   image(x,y,z,col=grey(seq(0.5,1,0.02)))
   contour(x,y,z,add=TRUE,col=1)
   mtext(side=3,line=2, namesvar[vcol]);
  }
  return(variab.ras)
}

make.variogram <- function (nugget=0, sill = 1000, range = 1000) {
    spherical.v <- function(h, parameters) ifelse(h == 0, 0, 
        ifelse(h <= parameters[3], parameters[1] + parameters[2] * 
            (3/2 * h/parameters[3] - 1/2 * (h/parameters[3])^3), 
            parameters[1] + parameters[2]))
    parameters <- c(nugget, sill, range)
    v.m.object <- list(parameters = parameters, model = spherical.v)
    names(v.m.object$parameters) <- c("nugget", "sill", "range")
    attr(v.m.object, "class") <- "variogram.model"
    attr(v.m.object, "type") <- "spherical"
    return(v.m.object)
}
 


