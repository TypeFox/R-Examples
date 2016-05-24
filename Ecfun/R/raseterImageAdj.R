rasterImageAdj <- function(image, xleft=par('usr')[1],
     ybottom=par('usr')[3], xright=par('usr')[2],
     ytop=par('usr')[4], angle = 0, interpolate = TRUE,
     xsub=NULL, ysub=NULL, ...){
##
## 1.  x, y pixels in image
##
  if(is.null(image)){
    stop("input 'image' is NULL")
  }
  Image <- as.raster(image)
#   dim(image) = y, x, RGB
  imagePixels <- dim(Image)[2:1]
#    imageAsp <- imagePixels[2]/imagePixels[1]
  names(imagePixels) <- c('x', 'y')
  if(!(is.null(xsub) && is.null(ysub))){
    if(is.null(xsub)){
      xsub <- 1:imagePixels[1]
    }
    if(is.null(ysub)){
      ysub <- 1:imagePixels[2]
    }
    if(length(dim(Image))>2){
      Image <- Image[ysub, xsub, , drop=FALSE]
    } else {
      Image <- Image[ysub, xsub, drop=FALSE]
    }    
    imagePixels[1:2] <- dim(Image)[2:1]
  }
##
## 2.  x, y pixels per inch
##
#  2.1.  x, y units in specified region
#   Can't compuute with x and y together:  
#   if one is a Date, it sometimes generates
#   a strange error below 
#    imageUnits <- c(x=xright-xleft, y=ytop-ybottom)
    imageUnits.x <- (xright-xleft)
    imageUnits.y <- (ytop-ybottom)
#    if(par('xlog'))imageUnits[1] <- log10(xright/xleft)
#    if(par('ylog'))imageUnits[2] <- log10(ytop/ybottom)
    if(par('xlog'))imageUnits.x <- log10(xright/xleft)
    if(par('ylog'))imageUnits.y <- log10(ytop/ybottom)
#  2.2.   plot units per inch
    xyinches <- xyinch(warn.log=FALSE)
#    plotAsp <- xyinches[2]/xyinches[1]
    names(xyinches) <- c('x', 'y')
#  2.3.  x, y pixels per inch in image region
#    pixelsPerInch <- (imagePixels*xyinches/imageUnits)
    pixelsPerInch.x <- (imagePixels[1]*xyinches[1] / 
                          as.numeric(imageUnits.x)) 
    pixelsPerInch.y <- (imagePixels[2]*xyinches[2] / 
                          imageUnits.y)
##
## 3.  Shrink imageUnits to max(PixelsPerInch)
##
#    imageUnitsAdj <- (imagePixels*xyinches/max(pixelsPerInch))
    maxPPI <- max(pixelsPerInch.x, pixelsPerInch.y)
    imageUAdj.x <- (imagePixels[1] * xyinches[1] / maxPPI) 
    imageUAdj.y <- (imagePixels[2] * xyinches[2] / maxPPI)
##
## 4.  (dX, dY) = imageUnitsAdj/2 
##              = half of the (width, height) in plotting units.  
##
#  dXY <- imageUnitsAdj/2 
  dX <- imageUAdj.x/2 
  dY <- imageUAdj.y/2 
# lower left = 45 degrees from center 
# Therefore, max deviation of a corner from the center:  at 45 degrees
#  dX. <- dXY[1]*sqrt(2)
#  dY. <- dXY[2]*sqrt(2)
  dX. <- dX*sqrt(2) 
  dY. <- dY*sqrt(2) 
##
## 5.  cntr = (xleft, ybottom) + imageUnits/2   
##
#  cntr <- (c(xleft, ybottom) + (imageUnits/2)) 
  cntr.x <- (xleft+(imageUnits.x/2))
  cntr.y <- (ybottom+(imageUnits.y/2))
  if(par('xlog')){
#    cntr[1] <- (xleft*10^(imageUnits[1]/2)) 
    cntr.x <- (xleft*10^(imageUnits.x/2)) 
  } 
  if(par('ylog')){
#    cntr[2] <- (ybottom*10^(imageUnits[2]/2))
    cntr.y <- (ybottom*10^(imageUnits.y/2))
  } 
##
## 6.  (x, y) location of the nominal lower left corner
##     after rotation 
##
  adj.x <- sin((angle-45)*pi/180)*dX. 
  adj.y <- cos((angle-45)*pi/180)*dY.
  if(par('xlog')){
#    xleft0 <- cntr[1]*10^adj.x
#    xright0 <- xleft0*10^imageUnitsAdj[1] 
    xleft0 <- cntr.x*10^adj.x
    xright0 <- xleft0*10^imageUAdj.x 
  } else {
    xleft0 <- (cntr.x + adj.x)
    xright0 <- (xleft0 + imageUAdj.x)
  }
  if(par('ylog')){
    ybottom0 <- cntr.y*10^(-adj.y)
    ytop0 <- (ybottom0*10^imageUAdj.y)
  } else {
    ybottom0 <- (cntr.y - adj.y) 
    ytop0 <- (ybottom0 + imageUAdj.y) 
  }
##
## 7.  rasterImage(image, xleft0, ybottom0, 
##        xright0, ytop0, angle, interpolate, ...)
##
  rasterImage(image=Image, xleft=xleft0, ybottom=ybottom0,
                xright=xright0, ytop=ytop0, angle=angle,
                interpolate=interpolate, ...)
##
## 8.  Done 
## 
  invisible(c(xleft=as.numeric(xleft0), 
    ybottom=as.numeric(ybottom0), xright=as.numeric(xright0), 
    ytop=as.numeric(ytop0)))
}
