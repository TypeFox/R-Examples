########## R-function: get.pixel.info ##########

# For settiing up pixel information for an
# image plot.

# Last changed: 18 JAN 2005

get.pixel.info <- function(plotControl.obj)
{
   # Unpack relevant components of plotControl.obj

   bdry <- plotControl.obj$bdry

   image.xlim <- plotControl.obj$image.xlim
   image.ylim <- plotControl.obj$image.ylim
   image.grid.size <- plotControl.obj$image.grid.size
   
   # Obtain grids in each direction according to
   # whether or not a boundary is specified.
 
   if (is.null(bdry))   # Use inputs to specify boundary
   {
      x.grid <- seq(image.xlim[1],image.xlim[2],length=image.grid.size[1])
      y.grid <- seq(image.ylim[1],image.ylim[2],length=image.grid.size[2])

      x.len <- length(x.grid)
      y.len <- length(y.grid)

      num.pix <- x.len*y.len
      pixel.mat <- expand.grid(x.grid,y.grid)
      
      pixels.on <- 1:num.pix
   }

   if (!is.null(bdry))  # Use "bdry" to specify pixels in absence of
   {                    # range information

      if (is.null(image.xlim))
         x.grid <- seq(min(bdry[,1]),max(bdry[,1]),length=image.grid.size[1])

      if (!is.null(image.xlim))
         x.grid <- seq(image.xlim[1],image.xlim[2],length=image.grid.size[1])

      if (is.null(image.ylim))
         y.grid <- seq(min(bdry[,2]),max(bdry[,2]),length=image.grid.size[2])

      if (!is.null(image.ylim))
         y.grid <- seq(image.ylim[1],image.ylim[2],length=image.grid.size[2])

      x.len <- length(x.grid)
      y.len <- length(y.grid)

      num.pix <- x.len*y.len

      pixel.mat <- expand.grid(x.grid,y.grid)

      pixels.on <- pointsInPoly(pixel.mat,bdry) 
   }

   x.frame <- pixel.frame(x.grid)
   y.frame <- pixel.frame(y.grid)

   on.inds <- (1:num.pix)[pixels.on]
   off.inds <- (1:num.pix)[!pixels.on]

   newdata <- as.data.frame(pixel.mat[on.inds,])

   return(list(x.frame=x.frame,y.frame=y.frame,on.inds=on.inds,
          off.inds=off.inds,newdata=newdata))
 
}

######### End of get.pixel.info ##########












