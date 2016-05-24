triangulate<-function(outer.poly,holes) {
	if(is.vector(outer.poly)&&is.numeric(outer.poly)&&(length(outer.poly)==4)) {
		stopifnot((outer.poly[3]-outer.poly[1])>0)
		stopifnot((outer.poly[4]-outer.poly[2])>0)
		outer.poly<-list(x=c(outer.poly[1],outer.poly[1],outer.poly[3],outer.poly[3]),
		y=c(outer.poly[2],outer.poly[4],outer.poly[4],outer.poly[2]))
	}
	if(!missing(holes)) {
		if(is.list(holes[[1]])) {
			if(all(unlist(lapply(holes,is.poly))))
				stopifnot(!overlapping.polygons(holes))
		}
		else if(is.poly(holes))
			holes<-list(holes)
		for(i in 1:length(holes))
			stopifnot(in.poly(holes[[i]]$x,holes[[i]]$y,outer.poly))
		nbpoly<-length(holes)+1
		nbpts<-length(outer.poly$x)
		vertX<-outer.poly$x
		vertY<-outer.poly$y
		for(i in 2:nbpoly) {
			nbpts[i]<-length(holes[[i-1]]$x)
			vertX<-c(vertX,holes[[i-1]]$x)
			vertY<-c(vertY,holes[[i-1]]$y)
		}
		nbptot<-sum(nbpts)
		nbtri<-(nbptot-2)+2*(nbpoly-1)
	}
	else {
		nbpoly<-1
		nbpts<-nbptot<-length(outer.poly$x)
		vertX<-outer.poly$x
		vertY<-outer.poly$y
		nbtri<-(nbpts-2)
	}
	tri<-.C("triangulate",
		as.integer(nbpoly),as.integer(nbpts),as.integer(nbptot),as.double(vertX),as.double(vertY),as.integer(nbtri),
		X1=double(nbtri),Y1=double(nbtri),X2=double(nbtri),Y2=double(nbtri),X3=double(nbtri),Y3=double(nbtri),
		PACKAGE="ads")
	return(data.frame(ax=tri$X1,ay=tri$Y1,bx=tri$X2,by=tri$Y2,cx=tri$X3,cy=tri$Y3))
}

