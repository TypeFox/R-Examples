dval<-function(p,upto,by,nx,ny) {
#si multivariŽ, choix du type de points ???
	stopifnot(inherits(p,"spp"))
	stopifnot(is.numeric(upto))
	stopifnot(is.numeric(by))
	stopifnot(by>0)
	r<-seq(by,upto,by)
	tmax<-length(r)
	stopifnot(is.numeric(nx))
	stopifnot(nx>=1)
	nx<-testInteger(nx)
	stopifnot(is.numeric(ny))
	stopifnot(ny>=1)
	ny<-testInteger(ny)		
	if("rectangle"%in%p$window$type) {
		cas<-1
		xmin<-p$window$xmin
		xmax<-p$window$xmax
		ymin<-p$window$ymin
		ymax<-p$window$ymax
		stopifnot(upto<=(0.5*max((xmax-xmin),(ymax-ymin))))
		xsample<-rep(xmin+(seq(1,nx)-0.5)*(xmax-xmin)/nx,each=ny)
		ysample<-rep(ymin+(seq(1,ny)-0.5)*(ymax-ymin)/ny,nx)
		if ("complex"%in%p$window$type) {
			cas<-3
			tri<-p$window$triangles
			nbTri<-nrow(tri)
		}
	}
	else if("circle"%in%p$window$type) {
		cas<-2
		x0<-p$window$x0
		y0<-p$window$y0
		r0<-p$window$r0
		stopifnot(upto<=r0)
		xsample<-rep(x0-r0+(seq(1,nx)-0.5)*2*r0/nx,each=ny)
		ysample<-rep(y0-r0+(seq(1,ny)-0.5)*2*r0/ny,nx)
		if ("complex"%in%p$window$type) {
			cas<-4
			tri<-p$window$triangles
			nbTri<-nrow(tri)
		}
	}
	else
		stop("invalid window type")
	
	ok <- inside.swin(xsample, ysample, p$window)
	xsample<-xsample[ok]
	ysample<-ysample[ok]
	stopifnot(length(xsample)==length(ysample))
	nbSample<-length(xsample)
	
	if(cas==1) { #rectangle
		count<-.C("density_rect",
				as.integer(p$n),as.double(p$x),as.double(p$y),
				as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
				as.integer(tmax),as.double(by),
				as.double(xsample),as.double(ysample),as.integer(nbSample),
				count=double(tmax*nbSample),
				PACKAGE="ads")$count
	}
	else if(cas==2) { #circle
		count<-.C("density_disq",
				as.integer(p$n),as.double(p$x),as.double(p$y),
				as.double(x0),as.double(y0),as.double(r0),
				as.integer(tmax),as.double(by),
				as.double(xsample),as.double(ysample),as.integer(nbSample),
				count=double(tmax*nbSample),
				PACKAGE="ads")$count
	}
	else if(cas==3) { #complex within rectangle
		count<-.C("density_tr_rect",
				as.integer(p$n),as.double(p$x),as.double(p$y),
				as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
				as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
				as.integer(tmax),as.double(by),
				as.double(xsample),as.double(ysample),as.integer(nbSample),
				count=double(tmax*nbSample),
				PACKAGE="ads")$count
	}
	else if(cas==4) { #complex within circle
		count<-.C("density_tr_disq",
				as.integer(p$n),as.double(p$x),as.double(p$y),
				as.double(x0),as.double(y0),as.double(r0),
				as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
				as.integer(tmax),as.double(by),
				as.double(xsample),as.double(ysample),as.integer(nbSample),
				count=double(tmax*nbSample),
				PACKAGE="ads")$count
	}
	## rajouter un indice lorsque les disques ne sont pas indŽpendants
	# formatting results	
	dens<-count/(pi*r^2)
	#grid<-matrix(c(xsample,ysample),nrow=nbSample,ncol=2)
	count<-matrix(count,nrow=nbSample,ncol=tmax,byrow=TRUE)
	dens<-matrix(dens,nrow=nbSample,ncol=tmax,byrow=TRUE)
	call<-match.call()	
	res<-list(call=call,window=p$window,r=r,xy=data.frame(x=xsample,y=ysample),cval=count,dval=dens)
	class(res)<-c("vads","dval")
	return(res)	
}

kval<-function(p,upto,by) {	
	# checking for input parameters
	stopifnot(inherits(p,"spp"))
	if(p$type!="univariate")
		warning(paste(p$type,"point pattern has been considered to be univariate\n"))
	stopifnot(is.numeric(upto))
	stopifnot(is.numeric(by))
	stopifnot(by>0)
	r<-seq(by,upto,by)
	tmax<-length(r)
	
	if("rectangle"%in%p$window$type) {
		cas<-1
		xmin<-p$window$xmin
		xmax<-p$window$xmax
		ymin<-p$window$ymin
		ymax<-p$window$ymax
		stopifnot(upto<=(0.5*max((xmax-xmin),(ymax-ymin))))
		if ("complex"%in%p$window$type) {
			cas<-3
			tri<-p$window$triangles
			nbTri<-nrow(tri)
		}
	}
	else if("circle"%in%p$window$type) {
		cas<-2
		x0<-p$window$x0
		y0<-p$window$y0
		r0<-p$window$r0
		stopifnot(upto<=r0)
		if ("complex"%in%p$window$type) {
			cas<-4
			tri<-p$window$triangles
			nbTri<-nrow(tri)
		}
	}
	else
		stop("invalid window type")
	intensity<-p$n/area.swin(p$window)

	#computing ripley local functions
	if(cas==1) { #rectangle
		res<-.C("ripleylocal_rect",		
				as.integer(p$n),as.double(p$x),as.double(p$y),
				as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
				as.integer(tmax),as.double(by),
				gi=double(p$n*tmax),ki=double(p$n*tmax),
				PACKAGE="ads")
	}
	else if(cas==2) { #circle
		res<-.C("ripleylocal_disq",
				as.integer(p$n),as.double(p$x),as.double(p$y),
				as.double(x0),as.double(y0),as.double(r0),
				as.integer(tmax),as.double(by),
				gi=double(p$n*tmax),ki=double(p$n*tmax),
				PACKAGE="ads")
	}
	else if(cas==3) { #complex within rectangle
		res<-.C("ripleylocal_tr_rect",
				as.integer(p$n),as.double(p$x),as.double(p$y),
				as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
				as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
				as.integer(tmax),as.double(by),
				gi=double(p$n*tmax),ki=double(p$n*tmax),
				PACKAGE="ads")
	}
	else if(cas==4) { #complex within circle
		res<-.C("ripleylocal_tr_disq",
				as.integer(p$n),as.double(p$x),as.double(p$y),
				as.double(x0),as.double(y0),as.double(r0),
				as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
				as.integer(tmax),as.double(by),
				gi=double(p$n*tmax),ki=double(p$n*tmax),
				PACKAGE="ads")
	}	
	# formatting results
	#coord<-matrix(c(X$x,X$y),nrow=nbPts,ncol=2)
	#coord<-data.frame(x=p$x,y=p$y)
	#r<-seq(dr,dr*tmax,dr)
	#ds<-pi*r^2-pi*seq(0,dr*tmax-dr,dr)^2
	ds<-c(pi,diff(pi*r^2))
	gi<-matrix(res$gi/(intensity*ds),nrow=p$n,ncol=tmax,byrow=TRUE)
	ni<-matrix(res$ki/(pi*r^2),nrow=p$n,ncol=tmax,byrow=TRUE)
	ki<-matrix(res$ki/intensity,nrow=p$n,ncol=tmax,byrow=TRUE)
	li<-matrix(sqrt(res$ki/(intensity*pi))-r,nrow=p$n,ncol=tmax,byrow=TRUE)
	call<-match.call()
	res<-list(call=call,window=p$window,r=r,xy=data.frame(x=p$x,y=p$y),gval=gi,kval=ki,nval=ni,lval=li)
	class(res)<-c("vads","kval")
	return(res)
}

k12val<-function(p,upto,by,marks) {
	# checking for input parameters
	stopifnot(inherits(p,"spp"))
	stopifnot(p$type=="multivariate")
	stopifnot(is.numeric(upto))
	stopifnot(is.numeric(by))
	stopifnot(by>0)
	r<-seq(by,upto,by)
	tmax<-length(r)
	if(missing(marks))
		marks<-c(1,2)
	stopifnot(length(marks)==2)
	stopifnot(marks[1]!=marks[2])
	mark1<-marks[1]
	mark2<-marks[2]
	if(is.numeric(mark1))
		mark1<-levels(p$marks)[testInteger(mark1)]
	else if(!mark1%in%p$marks) stop(paste("mark \'",mark1,"\' doesn\'t exist",sep=""))
	if(is.numeric(mark2))
		mark2<-levels(p$marks)[testInteger(mark2)]
	else if(!mark2%in%p$marks) stop(paste("mark \'",mark2,"\' doesn\'t exist",sep=""))
	# initializing variables
	if("rectangle"%in%p$window$type) {
		cas<-1
		xmin<-p$window$xmin
		xmax<-p$window$xmax
		ymin<-p$window$ymin
		ymax<-p$window$ymax
		stopifnot(upto<=(0.5*max((xmax-xmin),(ymax-ymin))))
		if ("complex"%in%p$window$type) {
			cas<-3
			tri<-p$window$triangles
			nbTri<-nrow(tri)
		}
	}
	else if("circle"%in%p$window$type) {
		cas<-2
		x0<-p$window$x0
		y0<-p$window$y0
		r0<-p$window$r0
		stopifnot(upto<=r0)
		if ("complex"%in%p$window$type) {
			cas<-4
			tri<-p$window$triangles
			nbTri<-nrow(tri)
		}
	}
	else
		stop("invalid window type")
	surface<-area.swin(p$window)
	x1<-p$x[p$marks==mark1]
	y1<-p$y[p$marks==mark1]
	x2<-p$x[p$marks==mark2]
	y2<-p$y[p$marks==mark2]
	nbPts1<-length(x1)
	nbPts2<-length(x2)
	intensity2<-nbPts2/surface
	#computing intertype local functions
	if(cas==1) { #rectangle
		res<-.C("intertypelocal_rect",		
				as.integer(nbPts1),as.double(x1),as.double(y1),
				as.integer(nbPts2),as.double(x2),as.double(y2),
				as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
				as.integer(tmax),as.double(by),
				gi=double(nbPts1*tmax),ki=double(nbPts1*tmax),
				PACKAGE="ads")
	}
	else if(cas==2) { #circle
		res<-.C("intertypelocal_disq",
				as.integer(nbPts1),as.double(x1),as.double(y1),
				as.integer(nbPts2),as.double(x2),as.double(y2),
				as.double(x0),as.double(y0),as.double(r0),
				as.integer(tmax),as.double(by),
				gi=double(nbPts1*tmax),ki=double(nbPts1*tmax),
				PACKAGE="ads")
	}
	else if(cas==3) { #complex within rectangle
		res<-.C("intertypelocal_tr_rect",		
				as.integer(nbPts1),as.double(x1),as.double(y1),
				as.integer(nbPts2),as.double(x2),as.double(y2),
				as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
				as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
				as.integer(tmax),as.double(by),
				gi=double(nbPts1*tmax),ki=double(nbPts1*tmax),
				PACKAGE="ads")
	}
	else if(cas==4) { #complex within circle
		res<-.C("intertypelocal_tr_disq",
				as.integer(nbPts1),as.double(x1),as.double(y1),
				as.integer(nbPts2),as.double(x2),as.double(y2),
				as.double(x0),as.double(y0),as.double(r0),
				as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
				as.integer(tmax),as.double(by),
				gi=double(nbPts1*tmax),ki=double(nbPts1*tmax),
				PACKAGE="ads")
	}
	# formatting results
	#coord<-matrix(c(x1,y1),nrow=nbPts1,ncol=2)
	#coord<-data.frame(x1=x1,y1=y1)
	#r<-seq(dr,dr*tmax,dr)
	#ds<-pi*r^2-pi*seq(0,dr*tmax-dr,dr)^2
	ds<-c(pi,diff(pi*r^2))
	gi<-matrix(res$gi/(intensity2*ds),nrow=nbPts1,ncol=tmax,byrow=TRUE)
	ni<-matrix(res$ki/(pi*r^2),nrow=nbPts1,ncol=tmax,byrow=TRUE)
	ki<-matrix(res$ki/intensity2,nrow=nbPts1,ncol=tmax,byrow=TRUE)
	li<-matrix(sqrt(res$ki/(intensity2*pi))-r,nrow=nbPts1,ncol=tmax,byrow=TRUE)
	call<-match.call()
	res<-list(call=call,window=p$window,r=r,xy=data.frame(x=x1,y=y1),g12val=gi,k12val=ki,n12val=ni,l12val=li,marks=c(mark1,mark2))
	class(res)<-c("vads","k12val")
	return(res)
}

