circsizer.map<-function(circsizer.object,type,zero,clockwise,title=NULL,
labels=NULL,label.pos=NULL,rad.pos=NULL,raw.data=FALSE){

	if(class(circsizer.object)!="circsizer") stop("Argument 'circsizer.object' must be the output of the
	'circsizer.density' or 'cirsizer.regression' functions")
	if (is.null(labels)){
		if (!any(type==1:5)){ 
			type=3
			warning("Value specified for argument 'type' is not valid. 'type=3' was used")
		}
		label.pos <- rad.pos <- seq(0,7*pi/4,by=pi/4)
	}else{ 
		if(length(labels)!=length(label.pos)){
			warning("Arguments 'labels' and 'label.pos' are not the same length. 'type=3' was used")
			type=3
			label.pos <- rad.pos <- seq(0,7*pi/4,by=pi/4)
		}else{
			type=0
			if (is.null(rad.pos)) rad.pos <- seq(0,7*pi/4,by=pi/4)
		}
	}
	if (type==5){label.pos <- seq(pi/12,23*pi/12,by=pi/6)
			 rad.pos <- seq(0,11*pi/6,by=pi/6) 
	}
	if (type==1) labels <- c("N","NE","E","SE","S","SW","W","NW")
	if (type==2) labels <- c("0h","3h","6h","9h","12h","15h","18h","21h")
	if (type==3) labels <- c(expression(0),expression(pi/4),expression(pi/2),expression(3*pi/4),expression(pi),
		                   expression(5*pi/4),expression(3*pi/2),expression(7*pi/4))
	if (type==4) labels <- c("0","45","90","135","180","225","270","315")
	if (type==5) labels <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Ago","Sep","Oct","Nov","Dec")
	
	# Create the plot
	bw.plot<-circsizer.object$bw
	nbw <- length(bw.plot) 
	stbw <- (bw.plot[2]-bw.plot[1])/2 
	radlim<- c((bw.plot[1]-stbw),bw.plot[nbw]+stbw)
	radlim[1]<-radlim[1]-diff(radlim)*0.2
	grid.pos<-pretty(radlim)
	st<-grid.pos[2]-grid.pos[1]
	radial.plot(0,rp.type="s",labels="",point.col="white",radial.lim=radlim,show.grid=F,show.radial.grid=F,
     	radial.labels="",show.grid.labels=F,start=zero,clockwise=clockwise,main=title) 
  
	# Add the labels
	label.prop=1.12
	maxlength=diff(radlim)
	if (clockwise) label.pos <- -label.pos + zero
	else label.pos <- label.pos + zero
	xpos <- cos(label.pos) * maxlength * label.prop
	ypos <- sin(label.pos) * maxlength * label.prop
	text(xpos, ypos, labels, cex = par("cex.axis"))

	# Plot the arrow
	arrow.pos <- seq(pi/4-0.13, pi/4+0.13, length=250)  
	xarrowpos <- cos(arrow.pos) * maxlength * 1.28
	yarrowpos <- sin(arrow.pos) * maxlength * 1.28
	points.default(xarrowpos,yarrowpos,type="l",lwd=2)
	if (clockwise) Arrows(xarrowpos[1], yarrowpos[1], xarrowpos[1]+(maxlength*0.02), yarrowpos[1]-(maxlength*0.02), type="simple", arr.length=0.5, segment=FALSE)
	else Arrows(xarrowpos[250], yarrowpos[250], xarrowpos[250]-(maxlength*0.02), yarrowpos[250]+(maxlength*0.02), type="simple", arr.length=0.5, segment=FALSE)

	# Add the color rings
	ngrid <- circsizer.object$ngrid
	t <- seq(0,2*pi,length=ngrid)
	rad<-seq(0,2*pi,length=250)
	for (i in 1:nbw){
		for (j in 1:ngrid){
			col<-circsizer.object$col[i,j]
			stt<-(t[2]-t[1])/2
			if (circsizer.object$log.scale) bwi<-bw.plot[i]-stbw
			else bwi<-ifelse(bw.plot[i]-stbw<0,0,bw.plot[i]-stbw)
			bws<-bw.plot[i]+stbw
			ti<-t[j]-stt
			ts<-t[j]+stt
			v1<-c(ti,bwi)
			v2<-c(ti,bws)
			v3<-c(ts,bwi)
			v4<-c(ts,bws)
			v1v3<-seq(v1[1],v3[1],by=0.01)
			v4v2<-sort(v1v3,decreasing=T)
			radial.plot(c(v1[2],rep(bwi,length(v1v3)),v3[2],v4[2],rep(bws,length(v4v2)),v2[2]),c(v1[1],v1v3,v3[1],v4[1],v4v2,v2[1]),rp.type="p",poly.col=col,line.col=col,radial.lim=radlim,add=T,start=zero,clockwise=clockwise)
		}
	}
	# Add the rings
	xcir <- grid.pos
	ncir <- length(xcir)
	ind<-c(which(xcir>bw.plot[nbw]+stbw),which(xcir<bw.plot[1]-stbw))
	if(length(ind)>0){
		xcir<- xcir[-ind]
		ncir<- ncir-length(ind)
	}
	for(i in 1:ncir){radial.plot(rep(xcir[i],250),rad,rp.type="p",radial.lim=radlim,add=T)} 
	# Add the radius
	segments(rep(0,9),rep(0,9),cos(rad.pos)*maxlength,sin(rad.pos)*maxlength) 
	# Add the values of the smoothing parameter along the radius
	xpos <- xcir - radlim[1]
	ypos <- rep(0, ncir)
	boxed.labels(xpos,ypos,as.character(xcir),border = FALSE,cex=0.7) 
	# Add the raw data
	if (raw.data){
		x <- circsizer.object$data
		attr(x, "class") <- attr(x, "circularp") <- NULL
		if (is.null(dim(x))){
			bins <- length(x)
			if (clockwise) x <- -x + zero
			else x <- x + zero
			x <- x%%(2 * pi)
			x[x >= 2 * pi] <- 2 * pi - 4 * .Machine$double.eps
			arc <- (2 * pi)/bins
			pos.bins <- (1 - 1/2) * arc - arc/2
			breaks <- seq(0, 2 * pi, length.out = (bins + 1))
			bins.count <- hist.default(x, breaks = breaks, plot = FALSE, right = TRUE)$counts
			mids <- seq(arc/2, 2 * pi - pi/bins, length = bins) + pos.bins[1]
			for (i in 1:bins) {
				if (bins.count[i] != 0) {
					for (j in 0:(bins.count[i] - 1)) {
						r <- 1 + j * 0.025
                  			z <- r * cos(mids[i])+1
                  			y <- r * sin(mids[i])
						points.default(z*maxlength-maxlength,y*maxlength, pch = 16, cex = 1)
					}
				}
			}
		}
	}
}