lines.regression.circular<-function(x, plot.type=c("circle", "line"), points.plot=FALSE, rp.type="p", type="l",
line.col=1, points.col="grey", points.pch=1, units=NULL, zero=NULL, clockwise=NULL, radial.lim=NULL, plot.info=NULL, ...){

	xcircularp <- attr(x$x, "circularp")
	ycircularp <- attr(x$y, "circularp")
    	if (is.null(xcircularp) && is.null(ycircularp)) 
 	stop("the component 'x' and/or the component 'y' of the object must be of class circular")
      plot.type <- match.arg(plot.type)
	if (is.circular(x$datax) && !is.circular(x$datay)){
		if (is.null(units)) units <- xcircularp$units
  		template <- xcircularp$template
		x$x <- conversion.circular(x$x, units = "radians", modulo = "2pi")
    		x$datax <- conversion.circular(x$datax, units = "radians", modulo = "2pi")
		attr(x$x, "class") <- attr(x$x, "circularp") <- NULL  
    		attr(x$datax, "class") <- attr(x$datax, "circularp") <- NULL

		if (plot.type=="line" & units == "degrees") {
           		x$x <- x$x/pi * 180
           		x$datax <- x$datax/pi * 180
        	}
        	if (plot.type=="line" & units == "hours") {
            	x$x <- x$x/pi * 12
            	x$datax <- x$datax/pi * 12
       	}

	} else if (is.circular(x$datax) && is.circular(x$datay)){
		if (plot.type=="circle") plot.type <- "torus"
  		template <- xcircularp$template
		if (is.null(units)) units <- xcircularp$units
		x$x <- conversion.circular(x$x, units = "radians", modulo = "2pi")
    		x$datax <- conversion.circular(x$datax, units = "radians", modulo = "2pi")
		attr(x$x, "class") <- attr(x$x, "circularp") <- NULL  
    		attr(x$datax, "class") <- attr(x$datax, "circularp") <- NULL

  		template <- ycircularp$template
		x$y <- conversion.circular(x$y, units = "radians", modulo = "2pi")
		x$datay <- conversion.circular(x$datay, units = "radians", modulo = "2pi")
		attr(x$y, "class") <- attr(x$y, "circularp") <- NULL  
		attr(x$datay, "class") <- attr(x$datay, "circularp") <- NULL

		x$datax[x$datax>pi]<-x$datax[x$datax>pi]-2*pi
		x$datay[x$datay>pi]<-x$datay[x$datay>pi]-2*pi
		x$x[x$x>pi]<-x$x[x$x>pi]-2*pi
		x$y[x$y>pi]<-x$y[x$y>pi]-2*pi

		if (plot.type=="line" & units == "degrees") {
           		x$x <- x$x/pi * 180
           		x$datax <- x$datax/pi * 180
           		x$y <- x$y/pi * 180
           		x$datay <- x$datay/pi * 180
        	}
        	if (plot.type=="line" & units == "hours") {
            	x$x <- x$x/pi * 12
            	x$datax <- x$datax/pi * 12
           		x$y <- x$y/pi * 12
           		x$datay <- x$datay/pi * 12
       	}

	} else if (!is.circular(x$datax) && is.circular(x$datay)){
		if (plot.type=="circle") plot.type <- "cylinder"
  		template <- ycircularp$template
		if (is.null(units)) units <- ycircularp$units
		x$y <- conversion.circular(x$y, units = "radians", modulo = "2pi")
		x$datay <- conversion.circular(x$datay, units = "radians", modulo = "2pi")		
		attr(x$y, "class") <- attr(x$y, "circularp") <- NULL  
		attr(x$datay, "class") <- attr(x$datay, "circularp") <- NULL

		if (plot.type=="line" & units == "degrees") {
           		x$y <- x$y/pi * 180
           		x$datay <- x$datay/pi * 180
        	}
       	if (plot.type=="line" & units == "hours") {
            	x$y <- x$y/pi * 12
            	x$datay <- x$datay/pi * 12
       	}
	}	
	if (plot.type == "line") {
       	xorder <- order(x$x)
       	x$x <- x$x[xorder]
        	x$y <- x$y[xorder]
        	lines.default(x$x, x$y, type = type, col=line.col, ...)
        	if (points.plot) points(x$datax, x$datay, col=points.col, pch=points.pch, ...)
	} else {
		if (plot.type=="torus"){
			xx <-cos(x$x)*(1+0.25*cos(x$y))
			yy <- sin(x$x)*(1+0.25*cos(x$y))
			zz <- 0.25*sin(x$y)
			lines3d(xx, yy, zz, col=line.col, ...)
			if (points.plot) {
				xx <- cos(x$datax)*(1+0.25*cos(x$datay))
				yy <- sin(x$datax)*(1+0.25*cos(x$datay))
				zz <- 0.25*sin(x$datay)
				points3d(xx, yy, zz, col=points.col)
			}
		} else if (plot.type=="cylinder"){
			R<- diff(range(x$datax))/8
			xx <- x$x
			yy <- R*cos(x$y)
			zz <- R*sin(x$y)
			lines3d(xx, yy, zz, col=line.col, ...)
			if (points.plot) {
				xx <- x$datax
				yy <- R*cos(x$datay)
				zz <- R*sin(x$datay)
				points3d(xx, yy, zz, col=points.col)
			}
		}else{
    			if (is.null(plot.info)) {
				if (is.null(radial.lim)) radial.lim <- range(c(x$datay,x$y))
        			if (is.null(zero)) {
            			if (template == "geographics" | template == "clock24") zero <- pi/2
            			else zero <- xcircularp$zero
        			}
        			if (is.null(clockwise)) {
            			if (template == "geographics" | template == "clock24") clockwise <- TRUE
            			else clockwise <- ifelse(xcircularp$rotation=="counter", FALSE, TRUE)
        			}
    			} else {
        			zero <- plot.info$zero
        			clockwise <- plot.info$clockwise
        			radial.lim <- plot.info$radial.lim
    			}
			radial.plot(x$y, x$x, rp.type=rp.type, line.col=line.col, start=zero, clockwise=clockwise, radial.lim=radial.lim, add=TRUE, ...) 
			if (points.plot) {
            		radial.plot(x$datay, x$datax, rp.type="s", start=zero, clockwise=clockwise, radial.lim=radial.lim, 
				point.col=points.col, point.symbols=points.pch, add=TRUE, ...)
			}
		}
	}
}
