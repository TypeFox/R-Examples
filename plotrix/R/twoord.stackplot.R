
twoord.stackplot <- function(lx, rx, ldata, rdata, lcol, rcol, ltype, rtype, 
	border, rylab, lylab, xlab, ..., incrylim=NULL,
	halfwidth=0.4, leftfront=FALSE, mar = c(5, 4, 4, 4))
{
	ltype <- sapply(ltype, function(x) match.arg(x, c("p","l","b","c","o","bar")))
	rtype <- sapply(rtype, function(x) match.arg(x, c("p","l","b","c","o","bar")))
	
	incrylim <- ifelse(is.null(incrylim), 0, as.numeric(incrylim))

	#convert to matrix	
	if(is.vector(ldata))
		ldata <- as.matrix(ldata)
	if(is.vector(rdata))
		rdata <- as.matrix(rdata)
		
	
	#some default parameters
	if(missing(border))
		border <- "grey80"
	if(missing(xlab))
		xlab <- "x"
	if(missing(rylab))
		rylab <- "right y values"		
	if(missing(lylab))
		lylab <- "left y values"	
	if(missing(lcol))
		lcol <- palette()[1:NCOL(ldata)]		
	if(missing(rcol))
		rcol <- palette()[1:NCOL(rdata)]		
		
	
	xlimits <- range(lx, rx)	
	
	#assume line plot for the left data, and barplot for the right data
    oldmar <- par("mar")
    par(mar = mar)
    	
	if(leftfront)
	{
		twoord.stackplot(rx, lx, rdata, ldata, rcol, lcol, rtype, ltype, 
			border, lylab, rylab, xlab, ..., incrylim=NULL,
			halfwidth=0.4, leftfront=FALSE, mar = c(5, 4, 4, 4))
		return(invisible())
	}
	
	#------------------------------------------------------------------------
	#left y-axis plot
	if(NCOL(ldata) > 1)
	{
	    lcol <- rep(lcol, length=NCOL(ldata))
	    ltype <- rep(ltype, length=NCOL(ldata))
	}    	
	if(any(ltype == "bar"))
	{	
		#------------------	
		lylimits <- range(ifelse(ldata <0, ldata, 0), rowSums(ldata))
		lylimits[1] <- ifelse(lylimits[1] > 0, lylimits[1]*(1-incrylim), lylimits[1]*(1+incrylim)) 
		lylimits[2] <- ifelse(lylimits[2] > 0, lylimits[2]*(1+incrylim), lylimits[2]*(1-incrylim)) 
		
		
		plot(0, type = "n", axes = FALSE, xlim = xlimits, ylim = lylimits, 
			ylab="", xlab=xlab, ...)
            
		xbottom <- par("usr")[1]
		xylim <- par("usr")
		ly <- ldata[,1]
		rect(lx-halfwidth, ifelse(ly < 0, ly, xbottom), lx+halfwidth, 
			ifelse(ly > 0, ly, 0), col=lcol[1], border=border, ...)
	
		if(NCOL(ldata) > 1)
			for(i in 2:NCOL(ldata))
			{
				ly <- ldata[,i]
				rect(lx-halfwidth, ifelse(ly < 0, ly, xbottom)+ldata[,i-1], 
					lx+halfwidth, ifelse(ly > 0, ly, 0)+ldata[,i-1], col=lcol[i], 
					border=border, ...)
			}	
		
	}else
	{
		#------------------
		lylimits <- range(ldata)
		lylimits[1] <- ifelse(lylimits[1] > 0, lylimits[1]*(1-incrylim), lylimits[1]*(1+incrylim)) 
		lylimits[2] <- ifelse(lylimits[2] > 0, lylimits[2]*(1+incrylim), lylimits[2]*(1-incrylim)) 
		
		plot(lx, ldata[, 1], xlim=xlimits, ylim=lylimits, col=lcol[1], 
			type=ltype[1], axes=FALSE, ylab="", xlab=xlab, ...)
		
		if(NCOL(ldata) > 1)
			for(i in 2:NCOL(ldata))
				lines(lx, ldata[, i], col=lcol[i], type=ltype[i], ...)
	}
	
	xylim <- par("usr")
	mtext(lylab, 2, 2, col = lcol[1])
		
	axis(1) #x axis
	axat <- axis(2, col=lcol[1], labels=FALSE) #left y axis
	abline(v=xylim[1], col=lcol[1])
	mtext(axat, 2, 1, at = axat, col = lcol[1])

	box()	
	

	#------------------------------------------------------------------------
	#right y-axis plot
	par(new=TRUE)
		
	if(NCOL(rdata) > 1)
	{
	    rcol <- rep(rcol, length=NCOL(rdata))
	    rtype <- rep(rtype, length=NCOL(rdata))
	}    	
	
	if(any(rtype == "bar"))
	{
		#------------------	
		rylimits <- range(ifelse(rdata <0, rdata, 0), rowSums(rdata))
		rylimits[1] <- ifelse(rylimits[1] > 0, rylimits[1]*(1-incrylim), rylimits[1]*(1+incrylim)) 
		rylimits[2] <- ifelse(rylimits[2] > 0, rylimits[2]*(1+incrylim), rylimits[2]*(1-incrylim)) 
		
		ry <- rdata[,1]
		plot(0, type = "n", axes = FALSE, xlim = xlimits, ylim = rylimits, 
			ylab="", xlab="", ...)
            
		xbottom <- par("usr")[1]
		xylim <- par("usr")
		rect(rx-halfwidth, ifelse(ry < 0, ry, xbottom), rx+halfwidth, 
			ifelse(ry > 0, ry, 0), col=rcol[1], border=border, ...)
	
		if(NCOL(rdata) > 1)
			for(i in 2:NCOL(rdata))
			{
				ry <- rdata[,i]
				rect(rx-halfwidth, ifelse(ry < 0, ry, xbottom)+rdata[,i-1], 
				rx+halfwidth, ifelse(ry > 0, ry, 0)+rdata[,i-1], col=rcol[i], 
				border=border, ...)

			}	
	}else
	{
		#------------------	
		rylimits <- range(rdata)
		rylimits[1] <- ifelse(rylimits[1] > 0, rylimits[1]*(1-incrylim), rylimits[1]*(1+incrylim)) 
		rylimits[2] <- ifelse(rylimits[2] > 0, rylimits[2]*(1+incrylim), rylimits[2]*(1-incrylim)) 
		
		plot(rx, rdata[, 1], xlim=xlimits, ylim=rylimits, col=rcol[1], 
			type=rtype[1], axes=FALSE, ylab="", xlab="", ...)
		
		if(NCOL(rdata) > 1)		
			for(i in 2:NCOL(rdata))
				lines(rx, rdata[, i], col=rcol[i], type=rtype[i], ...)
	}
	
	
	axat <- axis(4, col=rcol[1], labels=FALSE) #right y axis
	abline(v=xylim[1], col=rcol[1])
	mtext(axat, 4, 1, at = axat, col = rcol[1])

	mtext(rylab, 4, 2, col = rcol[1])
	
	par(mar = oldmar)   	
}

