# rectangle: c(xmin,ymin,xmax,ymax)
# circle: c(x0,y0,r0)
# triangles: mat(x1,y1,x2,y2,x3,y3)
swin<-function(window,triangles) {
	if(inherits(window,"swin")) {
		stopifnot("simple"%in%window$type)
		if(missing(triangles))
			return(window)
		else if("rectangle"%in%window$type)
			window<-c(window$xmin,window$ymin,window$xmax,window$ymax)
		else if("circle"%in%window$type)
			window<-c(window$x0,window$y0,window$r0)
		else
			stop("invalid window type")
	}
	stopifnot(is.numeric(window))
	stopifnot(length(window)%in%c(3,4))
	if(!missing(triangles)) {
		if(is.vector(triangles))
			triangles<-matrix(triangles,1,6)
		else
			triangles<-as.matrix(triangles)
		stopifnot(is.numeric(triangles))
		stopifnot(dim(triangles)[2]==6)
		dimnames(triangles)[[2]]<-c("ax","ay","bx","by","cx","cy")
		if(dim(triangles)[1]>1)
			stopifnot(!overlapping.polygons(convert(triangles)))
		triangles<-data.frame(triangles)
	}
	if(length(window)==4) {
		stopifnot((window[3]-window[1])>0)
		stopifnot((window[4]-window[2])>0)
		xmin<-window[1]
		ymin<-window[2]
		xmax<-window[3]
		ymax<-window[4]
		sw<-list(type=c("simple","rectangle"),xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax)  
		if(!missing(triangles)) {
			sw$type=c("complex","rectangle")
			stopifnot(unlist(lapply(convert(triangles),function(x) in.rectangle(x$x,x$y,xmin,ymin,xmax,ymax))))
			sw$triangles<-triangles
		}
	}
	else if(length(window)==3) {
		x0<-window[1]
		y0<-window[2]
		r0<-window[3]
		sw<-list(type=c("simple","circle"),x0=x0,y0=y0,r0=r0)		
		if(!missing(triangles)) {
			sw$type=c("complex","circle")
			stopifnot(unlist(lapply(convert(triangles),function(x) in.circle(x$x,x$y,x0,y0,r0))))
			sw$triangles<-triangles
		}
	}
	else
		stop("invalid input parameters")
	class(sw) <- "swin"
	return(sw)
}

print.swin<-function (x, ...) {
	cat("Sampling window:\n")
	str(x)
}

area.swin<-function (w) {
    stopifnot(inherits(w,"swin"))
	if("rectangle"%in%w$type)
		area<-(w$xmax-w$xmin)*(w$ymax-w$ymin)
	else if("circle"%in%w$type)
		area<-pi*w$r0^2
	else
		stop("invalid window type")
	if ("complex"%in%w$type) {
		tri<-w$triangles
		area.tri<-0
		for(i in 1:nrow(tri))
			area.tri<-area.tri+abs(area.poly(c(tri$ax[i],tri$bx[i],tri$cx[i]),c(tri$ay[i],tri$by[i],tri$cy[i])))
		area<-area-area.tri
	}
    return(area)
}

summary.swin<-function (object, ...) {
	res<-alist()
	res$type<-object$type
	if("rectangle"%in%object$type) {
		res$xrange<-c(object$xmin,object$xmax)
		res$yrange<-c(object$ymin,object$ymax)
	}
	else if("circle"%in%object$type) {
		res$xrange<-c(object$x0-object$r0,object$x0+object$r0)
		res$yrange<-c(object$y0-object$r0,object$y0+object$r0)
	}
	else
		stop("invalid window type")
	res$area<-area.swin(object)
	if("complex"%in%object$type) {
		res$nbtri<-nrow(object$triangles)
		if("rectangle"%in%object$type)
			res$area.init<-area.swin(swin(c(object$xmin,object$ymin,object$xmax,object$ymax)))
		else if("circle"%in%object$type)
			res$area.init<-area.swin(swin(c(object$x0,object$y0,object$r0)))
	}	
	class(res) <- "summary.swin"
    return(res)
}

print.summary.swin<-function (x,...) {
	cat(paste("Sampling window type:",x$type[1],x$type[2],"\n"))
	cat(paste("xrange: [",signif(x$xrange[1]),",",signif(x$xrange[2]),"]\n"))
	cat(paste("yrange: [",signif(x$yrange[1]),",",signif(x$yrange[2]),"]\n"))
	if("simple"%in%x$type)
		cat(paste("area:",signif(x$area),"\n"))
	else if("complex"%in%x$type) {
		cat(paste("initial",x$type[2],"area:",signif(x$area.init),"\n"))
		cat(paste("number of triangles removed:",x$nbtri,"\n"))
		cat(paste("actual complex window area:",signif(x$area),"\n"))
	}
	else
		stop("invalid window type")
}

plot.swin<-function (x,main,edge,scale=TRUE,add=FALSE,csize=1,...) {
	if(missing(main)) 
        main<-deparse(substitute(x))
	if(missing(edge))
		edge<-0
	#if(options()$device=="windows"&&sys.nframe()<=2)
	#	csize<-0.75*csize
	par(cex=csize)
	if("rectangle"%in%x$type) {
		rx<-c(x$xmin,x$xmax)
		ry<-c(x$ymin,x$ymax)
		if(edge>0) {
			rx<-c(rx[1]-edge,rx[2]+edge)
			ry<-c(ry[1]-edge,ry[2]+edge)
		}
		if(scale)
			plot(rx,ry,asp=1,main=main,type="n",axes=TRUE,frame.plot=FALSE,xlab="",ylab="",...)
		else
			plot(rx,ry,asp=1,main=main,type="n",axes=FALSE,xlab="",ylab="",...)
		polygon(c(x$xmin,x$xmin,x$xmax,x$xmax),c(x$ymin,x$ymax,x$ymax,x$ymin))
	}
	else if("circle"%in%x$type) {
		rx<-c(x$x0-x$r0,x$x0+x$r0)
		ry<-c(x$y0-x$r0,x$y0+x$r0)
		if(edge>0) {
			rx<-c(rx[1]-edge,rx[2]+edge)
			ry<-c(ry[1]-edge,ry[2]+edge)
		}
		if(scale)
			plot(rx,ry,asp=1,main=main,type="n",axes=TRUE,frame.plot=FALSE,xlab="",ylab="",...)
		else
			plot(rx,ry,asp=1,main=main,type="n",axes=FALSE,xlab="",ylab="",...)
		symbols(x$x0,x$y0,circles=x$r0,inches=FALSE,add=TRUE)
	}
	else
		stop("invalid window type")
	if("complex"%in%x$type)	{
		tri<-x$triangles
		for(i in 1:length(tri$ax)) {
			xi<-c(tri$ax[i],tri$bx[i],tri$cx[i])
			yi<-c(tri$ay[i],tri$by[i],tri$cy[i])
			polygon(xi,yi,col="grey",...)
			text(mean(xi),mean(yi),labels=as.character(i),cex=1)
		}
	}
}

inside.swin<-function(x,y,w,bdry=TRUE) {
	stopifnot(inherits(w,"swin"))
	stopifnot(length(x)==length(y))
	if("rectangle"%in%w$type)
		inside<-in.rectangle(x,y,w$xmin,w$ymin,w$xmax,w$ymax,bdry)
	else if("circle"%in%w$type)
		inside<-in.circle(x,y,w$x0,w$y0,w$r0,bdry)
	else
		stop("invalid window type")
	if("complex"%in%w$type) {
		tri<-w$triangles
		for(i in 1:nrow(tri)) 
			inside[in.triangle(x,y,tri$ax[i],tri$ay[i],tri$bx[i],tri$by[i],tri$cx[i],tri$cy[i])]<-FALSE
	}   
	return(inside)
}

owin2swin<-function(w) {
	stopifnot(inherits(w,"owin"))
	if(identical(w$type,c("rectangle")))
		sw<-swin(c(w$xrange[1],w$yrange[1],w$xrange[2],w$yrange[2]))
	else if(identical(w$type,c("polygonal"))) {
		if(length(w$bdry)==1) { #single polygon
			stopifnot(w$bdry[[1]]$hole==FALSE)
			wx<-border(w,0.1,outside=TRUE)
			outer.poly<-data.frame(x=c(rep(wx$xrange[1],2),rep(wx$xrange[2],2)),y=c(wx$yrange,wx$yrange[2:1]))
			tri<-triangulate(outer.poly,data.frame(w$bdry[[1]][1:2]))
			sw<-swin(c(wx$xrange[1],wx$yrange[1],wx$xrange[2],wx$yrange[2]),triangles=tri)
		}
		else { #polygon with holes
			stopifnot(w$bdry[[1]]$hole==FALSE)
			bb<-bounding.box.xy(w$bdry[[1]][1:2])
			if((bb$xrange==w$xrange)&&(bb$yrange==w$yrange)&&(area.owin(bb)==w$bdry[[1]]$area)) {	#first poly is rectangular window frame
				outer.poly<-data.frame(x=c(rep(w$xrange[1],2),rep(w$xrange[2],2)),y=c(w$yrange,w$yrange[2:1]))
				for(i in 2:length(w$bdry)) {
					stopifnot(w$bdry[[i]]$hole==TRUE)
					if(i==2)
						tri<-triangulate(w$bdry[[i]][1:2])
					else
						tri<-rbind(tri,triangulate(w$bdry[[i]][1:2]))
				}
				sw<-swin(c(w$xrange[1],w$yrange[1],w$xrange[2],w$yrange[2]),triangles=tri)
			}
			else { #first poly is a polygonal frame
				wx<-border(w,0.1,outside=TRUE)
				outer.poly<-data.frame(x=c(rep(wx$xrange[1],2),rep(wx$xrange[2],2)),y=c(wx$yrange,wx$yrange[2:1]))
				tri<-triangulate(outer.poly,w$bdry[[1]][1:2])
				for(i in 2:length(w$bdry)) {
					stopifnot(w$bdry[[i]]$hole==TRUE)
					tri<-rbind(tri,triangulate(w$bdry[[i]][1:2]))
				}
				sw<-swin(c(wx$xrange[1],wx$yrange[1],wx$xrange[2],wx$yrange[2]),triangles=tri)
			}
		}
	}
	else
	stop("non convertible 'owin' object")
	return(sw)
}


