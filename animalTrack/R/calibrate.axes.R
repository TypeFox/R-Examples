calibrate.axes <-
function(x, y, scale=1, xcal=NULL, ycal=NULL){
 
  # Use x-axis as the one to calibrate to
  xmax <- max(x) ; xmin <- min(x)
  ymax <- max(y) ; ymin <- min(y)
  
  xsf <- max( ((ymax-ymin)/(xmax-xmin)), 1)
  ysf <- max( ((xmax-xmin)/(ymax-ymin)), 1)
  
  xoff <- ( (xmax+xmin)/2 ) * xsf
  yoff <- ( (ymax+ymin)/2 ) * xsf
  
  xh <- x*xsf - xoff ; xhscaled <- (xh/ (diff(range(x))/2) ) * scale
  yh <- y*ysf - yoff ; yhscaled <- (yh/ (diff(range(y))/2) ) * scale
  
  if (is.null(xcal) == FALSE){xcal <- (xcal/ (diff(range(x))/2) ) * scale}
  if (is.null(ycal) == FALSE){ycal <- (ycal/ (diff(range(x))/2) ) * scale}
  
  return(list("xsf" = xsf, "ysf" = ysf, "xoff" = xoff, "yoff" = yoff, "xh" = xh, "yh"= yh, 
              "xhscaled" = xhscaled, "yhscaled" = yhscaled, "xcal" = xcal, "ycal" = ycal))
}

