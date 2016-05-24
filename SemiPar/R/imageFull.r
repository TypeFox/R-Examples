########## R-function: imageFull ##########

# For drawing an image plot with full
# accompaniments.

# Last changed: 18 JAN 2005

imageFull <- function(z,plotControl.obj)
{
   # Unpack relevant components of plotControl.obj

   image.col <- plotControl.obj$image.col
   image.bg <-  plotControl.obj$image.bg
   image.zlim <- plotControl.obj$image.zlim
   image.bty <- plotControl.obj$image.bty
   image.main <- plotControl.obj$image.main
   image.xlab <- plotControl.obj$image.xlab
   image.ylab <- plotControl.obj$image.ylab
   add.legend <- plotControl.obj$add.legend
   leg.loc <- plotControl.obj$leg.loc
   leg.dim <- plotControl.obj$leg.dim 
   image.zlab.col <- plotControl.obj$image.zlab.col
   image.zlab <- plotControl.obj$image.zlab
 
   # Obtain relevant pixel information

   pixel.info <- get.pixel.info(plotControl.obj)

   x.frame <- pixel.info$x.frame
   y.frame <- pixel.info$y.frame

   # Determine default value of "image.zlim" parameter.

   clean.z <- z[is.na(z)==0]

   if (is.null(image.zlim)) 
       image.zlim <- range(clean.z)

   # Truncate values of z to lie within image.zlim

   z[z<image.zlim[1]] <- image.zlim[1]
   z[z>image.zlim[2]] <- image.zlim[2]
   
   # Draw image

   op <- par(bg="white")
   par(bg=image.bg)
   image(x.frame,y.frame,z,zlim=image.zlim,bty=image.bty,
         main=image.main,xlab=image.xlab,ylab=image.ylab,
         col=image.col)
   par(op)

   if (add.legend==TRUE)
   {
      
      # Add a legend to the plot, starting with the
      # formulation of the default location (if required)

      if (is.null(leg.loc))
      {
         leg.loc[1] <- min(x.frame) + 0.70*(max(x.frame)-min(x.frame))
         leg.loc[2] <- min(y.frame) + 0.10*(max(y.frame)-min(y.frame))
      }

      text.loc.x <- leg.loc[1] + 0.15*(max(x.frame)-min(x.frame))
      text.loc.y <- leg.loc[2] + 0.02*(max(x.frame)-min(x.frame))

      # Make sure that legend is not located too close to the
      # boundary of the image

      imageLegend(image.zlim,image.col,leg.loc,leg.dim,
                   image.zlab,image.zlab.col)
   }

   return(image.zlim)
}

######### End of imageFull ##########












