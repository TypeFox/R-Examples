kfun<-function(p,upto,by,nsim=0,prec=0.01,alpha=0.01) {
	# checking for input parameters
	stopifnot(inherits(p,"spp"))
	if(p$type!="univariate")
		warning(paste(p$type,"point pattern has been considered to be univariate\n"))
	stopifnot(is.numeric(upto))
	stopifnot(upto>=1)
	stopifnot(is.numeric(by))
	stopifnot(by>0)
	r<-seq(by,upto,by)
	tmax<-length(r)
	stopifnot(is.numeric(nsim))
	stopifnot(nsim>=0)
	nsim<-testInteger(nsim)	
	stopifnot(is.numeric(prec))
	stopifnot(prec>=0)
	stopifnot(is.numeric(alpha))
	stopifnot(alpha>=0)
	if(nsim>0) testIC(nsim,alpha)

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
	
	if(cas==1) { #rectangle
		if(nsim==0) { #without CI
			res<-.C("ripley_rect",
					as.integer(p$n),as.double(p$x),as.double(p$y),
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
					as.integer(tmax),as.double(by),
					g=double(tmax),k=double(tmax),
					PACKAGE="ads")
		}
		else { #with CI
			res<-.C("ripley_rect_ic",
					as.integer(p$n),as.double(p$x),as.double(p$y),
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),as.double(intensity),
					as.integer(tmax),as.double(by),
					as.integer(nsim),as.double(prec),as.double(alpha),
					g=double(tmax),k=double(tmax),
					gic1=double(tmax),gic2=double(tmax),kic1=double(tmax),kic2=double(tmax),
					gval=double(tmax),kval=double(tmax),lval=double(tmax),nval=double(tmax),
					PACKAGE="ads")
		}
	}
	else if(cas==2) { #circle
		if(nsim==0) { #without CI
			res<-.C("ripley_disq",
					as.integer(p$n),as.double(p$x),as.double(p$y),
					as.double(x0),as.double(y0),as.double(r0),
					as.integer(tmax),as.double(by),
					g=double(tmax),k=double(tmax),
					PACKAGE="ads")					
		}
		else { #with CI
			res<-.C("ripley_disq_ic",
					as.integer(p$n),as.double(p$x),as.double(p$y),
					as.double(x0),as.double(y0),as.double(r0),as.double(intensity),
					as.integer(tmax),as.double(by),
					as.integer(nsim),as.double(prec),as.double(alpha),
					g=double(tmax),k=double(tmax),
					gic1=double(tmax),gic2=double(tmax),kic1=double(tmax),kic2=double(tmax),
					gval=double(tmax),kval=double(tmax),lval=double(tmax),nval=double(tmax),
					PACKAGE="ads")
		}	
	}
	else if(cas==3) { #complex within rectangle
		if(nsim==0) { #without CI
			res<-.C("ripley_tr_rect",
					as.integer(p$n),as.double(p$x),as.double(p$y),
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by),
					g=double(tmax),k=double(tmax),
					PACKAGE="ads")					
		}
		else { #with CI
			res<-.C("ripley_tr_rect_ic",
					as.integer(p$n),as.double(p$x),as.double(p$y),
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),as.double(intensity),
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by),
					as.integer(nsim),as.double(prec),as.double(alpha),
					g=double(tmax),k=double(tmax),
					gic1=double(tmax),gic2=double(tmax),kic1=double(tmax),kic2=double(tmax),
					gval=double(tmax),kval=double(tmax),lval=double(tmax),nval=double(tmax),
					PACKAGE="ads")	
		}	
	}
	else if(cas==4) { #complex within circle
		if(nsim==0) { #without CI		
			res<-.C("ripley_tr_disq",
					as.integer(p$n),as.double(p$x),as.double(p$y),
					as.double(x0),as.double(y0),as.double(r0),
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by),
					g=double(tmax),k=double(tmax),
					PACKAGE="ads")	
		}
		else { #with CI
			res<-.C("ripley_tr_disq_ic",
					as.integer(p$n),as.double(p$x),as.double(p$y),
					as.double(x0),as.double(y0),as.double(r0),as.double(intensity),
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by),
					as.integer(nsim),as.double(prec),as.double(alpha),
					g=double(tmax),k=double(tmax),
					gic1=double(tmax),gic2=double(tmax),kic1=double(tmax),kic2=double(tmax),
					gval=double(tmax),kval=double(tmax),lval=double(tmax),nval=double(tmax),
					PACKAGE="ads")	
		}
	}
	# formatting results
	ds<-c(pi,diff(pi*r^2))
	g<-data.frame(obs=res$g/(intensity*ds),theo=rep(1,tmax))
	n<-data.frame(obs=res$k/(pi*r^2),theo=rep(intensity,tmax))
	k<-data.frame(obs=res$k/intensity,theo=pi*r^2)
	l<-data.frame(obs=sqrt(res$k/(intensity*pi))-r,theo=rep(0,tmax))
	if(nsim>0) {
		g<-cbind(g,sup=res$gic1/(intensity*ds),inf=res$gic2/(intensity*ds),pval=res$gval/(nsim+1))
		n<-cbind(n,sup=res$kic1/(pi*r^2),inf=res$kic2/(pi*r^2),pval=res$nval/(nsim+1))
		k<-cbind(k,sup=res$kic1/intensity,inf=res$kic2/intensity,pval=res$kval/(nsim+1))
		l<-cbind(l,sup=sqrt(res$kic1/(intensity*pi))-r,inf=sqrt(res$kic2/(intensity*pi))-r,pval=res$lval/(nsim+1))
	}
	call<-match.call()
	res<-list(call=call,r=r,g=g,n=n,k=k,l=l)
	class(res)<-c("fads","kfun")
	return(res)
}

k12fun<-function(p,upto,by,nsim=0,H0=c("pitor","pimim","rl"),prec=0.01,nsimax=3000,conv=50,rep=10,alpha=0.01,marks) {
	# checking for input parameters
	options( CBoundsCheck = TRUE )
	# regle les problemes pour 32-bit
	stopifnot(inherits(p,"spp"))
	stopifnot(p$type=="multivariate")
	stopifnot(is.numeric(upto))
	stopifnot(upto>=1)
	stopifnot(is.numeric(by))
	stopifnot(by>0)
	r<-seq(by,upto,by)
	tmax<-length(r)
	stopifnot(is.numeric(nsim))
	stopifnot(nsim>=0)
	nsim<-testInteger(nsim)
	H0<-H0[1]
	stopifnot(H0=="pitor" || H0=="pimim" || H0=="rl")
	if(H0=="rl") H0<-1
	else if(H0=="pitor") H0<-2
	else H0<-3
	stopifnot(is.numeric(prec))
	stopifnot(prec>=0)
	stopifnot(is.numeric(alpha))
	stopifnot(alpha>=0)
	if(nsim>0) testIC(nsim,alpha)
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
#	intensity<-(nbPts1+nbPts2)/surface
	
	# computing intertype functions
	if(cas==1) { #rectangle
		if(nsim==0) { #without CI
			res<-.C("intertype_rect",
					as.integer(nbPts1),as.double(x1),as.double(y1),
					as.integer(nbPts2),as.double(x2),as.double(y2),					
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),					
					as.integer(tmax),as.double(by),
					g=double(tmax),k=double(tmax),
					PACKAGE="ads")					
		}
		else { #with CI
			res<-.C("intertype_rect_ic",
					as.integer(nbPts1),as.double(x1),as.double(y1),
					as.integer(nbPts2),as.double(x2),as.double(y2),
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),as.double(surface),
					as.integer(tmax),as.double(by),
					as.integer(nsim),as.integer(H0),as.double(prec),as.integer(nsimax),as.integer(conv),as.integer(rep),as.double(alpha),
					g=double(tmax),k=double(tmax),
					gic1=double(tmax),gic2=double(tmax),kic1=double(tmax),kic2=double(tmax),
					gval=double(tmax),kval=double(tmax),lval=double(tmax),nval=double(tmax),
					PACKAGE="ads")
		}
	}
	else if(cas==2) { #circle
		if(nsim==0) { #without CI
			res<-.C("intertype_disq",
					as.integer(nbPts1),as.double(x1),as.double(y1),
					as.integer(nbPts2),as.double(x2),as.double(y2),
					as.double(x0),as.double(y0),as.double(r0),					
					as.integer(tmax),as.double(by),
					g=double(tmax),k=double(tmax),
					PACKAGE="ads")					
		}
		else { #with CI
			res<-.C("intertype_disq_ic",
					as.integer(nbPts1),as.double(x1),as.double(y1),
					as.integer(nbPts2),as.double(x2),as.double(y2),
					as.double(x0),as.double(y0),as.double(r0),as.double(surface),
					as.integer(tmax),as.double(by),
					as.integer(nsim),as.integer(H0),as.double(prec),as.integer(nsimax),as.integer(conv),as.integer(rep),as.double(alpha),
					g=double(tmax),k=double(tmax),
					gic1=double(tmax),gic2=double(tmax),kic1=double(tmax),kic2=double(tmax),
					gval=double(tmax),kval=double(tmax),lval=double(tmax),nval=double(tmax),
					PACKAGE="ads")		
		}
	}
	else if(cas==3) { #complex within rectangle
		if(nsim==0) { #without CI
			res<-.C("intertype_tr_rect",
					as.integer(nbPts1),as.double(x1),as.double(y1),
					as.integer(nbPts2),as.double(x2),as.double(y2),
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),					
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by),
					g=double(tmax),k=double(tmax),
					PACKAGE="ads")					
		}
		else { #with CI
			res<-.C("intertype_tr_rect_ic",
					as.integer(nbPts1),as.double(x1),as.double(y1),
					as.integer(nbPts2),as.double(x2),as.double(y2),
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),as.double(surface),
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by),
					as.integer(nsim),as.integer(H0),as.double(prec),as.integer(nsimax),as.integer(conv),as.integer(rep),as.double(alpha),
					g=double(tmax),k=double(tmax),
					gic1=double(tmax),gic2=double(tmax),kic1=double(tmax),kic2=double(tmax),
					gval=double(tmax),kval=double(tmax),lval=double(tmax),nval=double(tmax),
					PACKAGE="ads")		
		}
	}
	else if(cas==4) { #complex within circle
		if(nsim==0) { #without CI
			res<-.C("intertype_tr_disq",
					as.integer(nbPts1),as.double(x1),as.double(y1),
					as.integer(nbPts2),as.double(x2),as.double(y2),
					as.double(x0),as.double(y0),as.double(r0),					
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by),
					g=double(tmax),k=double(tmax),
					PACKAGE="ads")					
		}
		else { #with CI
			res<-.C("intertype_tr_disq_ic",
					as.integer(nbPts1),as.double(x1),as.double(y1),
					as.integer(nbPts2),as.double(x2),as.double(y2),
					as.double(x0),as.double(y0),as.double(r0),as.double(surface),
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by),
					as.integer(nsim),as.integer(H0),as.double(prec),as.integer(nsimax),as.integer(conv),as.integer(rep),as.double(alpha),
					g=double(tmax),k=double(tmax),
					gic1=double(tmax),gic2=double(tmax),kic1=double(tmax),kic2=double(tmax),
					gval=double(tmax),kval=double(tmax),lval=double(tmax),nval=double(tmax),
					PACKAGE="ads")		
		}
	}	
	# formatting results
	ds<-c(pi,diff(pi*r^2))
	g<-res$g/(intensity2*ds)
	n<-res$k/(pi*r^2)
	k<-res$k/intensity2
	l<-sqrt(res$k/(intensity2*pi))-r
	if(H0==1) {
		rip<-kfun(spp(c(x1,x2),c(y1,y2),p$window),upto,by)
		theo<-list(g=rip$g$obs,n=intensity2*rip$k$obs/(pi*r^2),k=rip$k$obs,l=rip$l$obs)
	}
	else if (H0==2||H0==3)
		theo<-list(g=rep(1,tmax),n=rep(intensity2,tmax),k=pi*r^2,l=rep(0,tmax))
	g<-data.frame(obs=g,theo=theo$g)
	n<-data.frame(obs=n,theo=theo$n)
	k<-data.frame(obs=k,theo=theo$k)
	l<-data.frame(obs=l,theo=theo$l)
	if(nsim>0) {
		g<-cbind(g,sup=res$gic1/(intensity2*ds),inf=res$gic2/(intensity2*ds),pval=res$gval/(nsim+1))
		n<-cbind(n,sup=res$kic1/(pi*r^2),inf=res$kic2/(pi*r^2),pval=res$nval/(nsim+1))
		k<-cbind(k,sup=res$kic1/intensity2,inf=res$kic2/intensity2,pval=res$kval/(nsim+1))
		l<-cbind(l,sup=sqrt(res$kic1/(intensity2*pi))-r,inf=sqrt(res$kic2/(intensity2*pi))-r,pval=res$lval/(nsim+1))
	}	
	call<-match.call()
	res<-list(call=call,r=r,g12=g,n12=n,k12=k,l12=l,marks=c(mark1,mark2))
	class(res)<-c("fads","k12fun")
	return(res)
}

kijfun<-kpqfun<-function(p,upto,by) {
# checking for input parameters
	stopifnot(inherits(p,"spp"))
	stopifnot(p$type=="multivariate")
	stopifnot(is.numeric(upto))
	stopifnot(upto>=1)
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
	surface<-area.swin(p$window)
	
	tabMarks<-levels(p$marks)
	nbMarks<-length(tabMarks)
	mpt_nb<-summary(p$marks)	
# computing RipleyClass
	gij<-double(tmax*nbMarks^2)
	kij<-double(tmax*nbMarks^2)
	lij<-double(tmax*nbMarks^2)
	nij<-double(tmax*nbMarks^2)
	for(i in 1:nbMarks) {
		x1<-p$x[p$marks==tabMarks[i]]
		y1<-p$y[p$marks==tabMarks[i]]		
		if(cas==1) { #rectangle
			res<-.C("ripley_rect",
					as.integer(mpt_nb[i]),as.double(x1),as.double(y1),
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
					as.integer(tmax),as.double(by),
					g=double(tmax),k=double(tmax),
					PACKAGE="ads")
		}
		else if(cas==2) { #circle
			res<-.C("ripley_disq",
					as.integer(mpt_nb[i]),as.double(x1),as.double(y1),
					as.double(x0),as.double(y0),as.double(r0),
					as.integer(tmax),as.double(by),
					g=double(tmax),k=double(tmax),
					PACKAGE="ads")
		}
		else if(cas==3) { #complex within rectangle
			res<-.C("ripley_tr_rect",
					as.integer(mpt_nb[i]),as.double(x1),as.double(y1),
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by),
					g=double(tmax),k=double(tmax),
					PACKAGE="ads")
		}
		else if(cas==4) { #complex within circle
			res<-.C("ripley_tr_disq",
					as.integer(mpt_nb[i]),as.double(x1),as.double(y1),
					as.double(x0),as.double(y0),as.double(r0),
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by),
					g=double(tmax),k=double(tmax),
					PACKAGE="ads")
		}
		intensity1<-mpt_nb[i]/surface
		matcol<-(i-1)*nbMarks+i-1
		j<-(matcol*tmax+1):(matcol*tmax+tmax)
		ds<-c(pi,diff(pi*r^2))
		gij[j]<-res$g/(intensity1*ds)
		nij[j]<-res$k/(pi*r^2)
		kij[j]<-res$k/intensity1
		lij[j]<-sqrt(res$k/(intensity1*pi))-r
		if(i<nbMarks) {
			for(j in (i+1):nbMarks) {
				x2<-p$x[p$marks==tabMarks[j]]
				y2<-p$y[p$marks==tabMarks[j]]		
				if(cas==1) { #rectangle
					res<-.C("intertype_rect",
							as.integer(mpt_nb[i]),as.double(x1),as.double(y1),
							as.integer(mpt_nb[j]),as.double(x2),as.double(y2),					
							as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),					
							as.integer(tmax),as.double(by),
							g=double(tmax),k=double(tmax),
							PACKAGE="ads")					
				}
				else if(cas==2) { #circle
					res<-.C("intertype_disq",
							as.integer(mpt_nb[i]),as.double(x1),as.double(y1),
							as.integer(mpt_nb[j]),as.double(x2),as.double(y2),
							as.double(x0),as.double(y0),as.double(r0),					
							as.integer(tmax),as.double(by),
							g=double(tmax),k=double(tmax),
							PACKAGE="ads")					
				}
				else if(cas==3) { #complex within rectangle
					res<-.C("intertype_tr_rect",
							as.integer(mpt_nb[i]),as.double(x1),as.double(y1),
							as.integer(mpt_nb[j]),as.double(x2),as.double(y2),
							as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),					
							as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
							as.integer(tmax),as.double(by),
							g=double(tmax),k=double(tmax),
							PACKAGE="ads")
				}
				else if(cas==4) { #complex within circle
					res<-.C("intertype_tr_disq",
							as.integer(mpt_nb[i]),as.double(x1),as.double(y1),
							as.integer(mpt_nb[j]),as.double(x2),as.double(y2),
							as.double(x0),as.double(y0),as.double(r0),					
							as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
							as.integer(tmax),as.double(by),
							g=double(tmax),k=double(tmax),
							PACKAGE="ads")					
				}
				intensity2<-mpt_nb[j]/surface
				matcol<-(i-1)*nbMarks+j-1
				k<-(matcol*tmax+1):(matcol*tmax+tmax)
				gij[k]<-res$g/(intensity2*ds)
				nij[k]<-res$k/(pi*r^2)
				kij[k]<-res$k/intensity2
				lij[k]<-sqrt(res$k/(intensity2*pi))-r				
				if(cas==1) { #rectangle
					res<-.C("intertype_rect",
							as.integer(mpt_nb[j]),as.double(x2),as.double(y2),
							as.integer(mpt_nb[i]),as.double(x1),as.double(y1),											
							as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),					
							as.integer(tmax),as.double(by),
							g=double(tmax),k=double(tmax),
							PACKAGE="ads")					
				}
				else if(cas==2) { #circle
					res<-.C("intertype_disq",
							as.integer(mpt_nb[j]),as.double(x2),as.double(y2),
							as.integer(mpt_nb[i]),as.double(x1),as.double(y1),						
							as.double(x0),as.double(y0),as.double(r0),					
							as.integer(tmax),as.double(by),
							g=double(tmax),k=double(tmax),
							PACKAGE="ads")					
				}
				else if(cas==3) { #complex within rectangle
					res<-.C("intertype_tr_rect",
							as.integer(mpt_nb[j]),as.double(x2),as.double(y2),
							as.integer(mpt_nb[i]),as.double(x1),as.double(y1),						
							as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),					
							as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
							as.integer(tmax),as.double(by),
							g=double(tmax),k=double(tmax),
							PACKAGE="ads")
				}
				else if(cas==4) { #complex within circle
					res<-.C("intertype_tr_disq",
							as.integer(mpt_nb[j]),as.double(x2),as.double(y2),
							as.integer(mpt_nb[i]),as.double(x1),as.double(y1),						
							as.double(x0),as.double(y0),as.double(r0),					
							as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
							as.integer(tmax),as.double(by),
							g=double(tmax),k=double(tmax),
							PACKAGE="ads")					
				}								
				matcol<-(j-1)*nbMarks+i-1
				k<-(matcol*tmax+1):(matcol*tmax+tmax)
				gij[k]<-res$g/(intensity1*ds)
				nij[k]<-res$k/(pi*r^2)
				kij[k]<-res$k/intensity1
				lij[k]<-sqrt(res$k/(intensity1*pi))-r
			}
		}
	}
	labij<-paste(rep(tabMarks,each=nbMarks),rep(tabMarks,nbMarks),sep="-")
	gij<-matrix(gij,nrow=tmax,ncol=nbMarks^2)
	kij<-matrix(kij,nrow=tmax,ncol=nbMarks^2)
	nij<-matrix(nij,nrow=tmax,ncol=nbMarks^2)
	lij<-matrix(lij,nrow=tmax,ncol=nbMarks^2)
	call<-match.call()
	res<-list(call=call,r=r,labpq=labij,gij=gij,kpq=kij,lpq=lij,npq=nij,intensity=summary(p)$intensity)
	class(res)<-c("fads","kpqfun")
	return(res)
}

ki.fun<-kp.fun<-function(p,upto,by) {
	# checking for input parameters
	stopifnot(inherits(p,"spp"))
		stopifnot(p$type=="multivariate")
	stopifnot(is.numeric(upto))
	stopifnot(upto>=1)
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
	surface<-area.swin(p$window)

	tabMarks<-levels(p$marks)
	nbMarks<-length(tabMarks)
	mpt_nb<-summary(p$marks)		
	#computing RipleyAll
	gis<-double(tmax*nbMarks)
	kis<-double(tmax*nbMarks)
	lis<-double(tmax*nbMarks)
	nis<-double(tmax*nbMarks)
	for(i in 1:nbMarks) {
		x1<-p$x[p$marks==tabMarks[i]]
		y1<-p$y[p$marks==tabMarks[i]]
		x2<-p$x[p$marks!=tabMarks[i]]
		y2<-p$y[p$marks!=tabMarks[i]]
		nbPts1<-mpt_nb[i]
		nbPts2<-sum(mpt_nb)-nbPts1		
		if(cas==1) { #rectangle
			res<-.C("intertype_rect",
					as.integer(nbPts1),as.double(x1),as.double(y1),
					as.integer(nbPts2),as.double(x2),as.double(y2),					
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),					
					as.integer(tmax),as.double(by),
					g=double(tmax),k=double(tmax),
					PACKAGE="ads")					
		}
		else if(cas==2) { #circle
			res<-.C("intertype_disq",
					as.integer(nbPts1),as.double(x1),as.double(y1),
					as.integer(nbPts2),as.double(x2),as.double(y2),
					as.double(x0),as.double(y0),as.double(r0),					
					as.integer(tmax),as.double(by),
					g=double(tmax),k=double(tmax),
					PACKAGE="ads")					
		}
		else if(cas==3) { #complex within rectangle
			res<-.C("intertype_tr_rect",
					as.integer(nbPts1),as.double(x1),as.double(y1),
					as.integer(nbPts2),as.double(x2),as.double(y2),
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),					
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by),
					g=double(tmax),k=double(tmax),
					PACKAGE="ads")
		}
		else if(cas==4) { #complex within circle
			res<-.C("intertype_tr_disq",
					as.integer(nbPts1),as.double(x1),as.double(y1),
					as.integer(nbPts2),as.double(x2),as.double(y2),
					as.double(x0),as.double(y0),as.double(r0),					
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by),
					g=double(tmax),k=double(tmax),
					PACKAGE="ads")					
		}		
		intensity2<-nbPts2/surface
		j<-((i-1)*tmax+1):((i-1)*tmax+tmax)
		ds<-c(pi,diff(pi*r^2))
		gis[j]<-res$g/(intensity2*ds)
		nis[j]<-res$k/(pi*r^2)
		kis[j]<-res$k/intensity2
		lis[j]<-sqrt(res$k/(intensity2*pi))-r
	}	
	# formatting results
	labi<-tabMarks
	gis<-matrix(gis,nrow=tmax,ncol=nbMarks)
	kis<-matrix(kis,nrow=tmax,ncol=nbMarks)
	nis<-matrix(nis,nrow=tmax,ncol=nbMarks)
	lis<-matrix(lis,nrow=tmax,ncol=nbMarks)
	call<-match.call()
	res<-list(call=call,r=r,labp=labi,gp.=gis,kp.=kis,lp.=lis,np.=nis,intensity=summary(p)$intensity)
	class(res)<-c("fads","kp.fun")
	return(res)
}

kmfun<-function(p,upto,by,nsim=0,alpha=0.01) {
	# checking for input parameters
	stopifnot(inherits(p,"spp"))
	stopifnot(p$type=="marked")
	stopifnot(is.numeric(upto))
	stopifnot(upto>=1)
	stopifnot(is.numeric(by))
	stopifnot(by>0)
	r<-seq(by,upto,by)
	tmax<-length(r)
	stopifnot(is.numeric(nsim))
	stopifnot(nsim>=0)
	nsim<-testInteger(nsim)
	#cmoy<-mean(p$marks)
	cvar<-var(p$marks)
	stopifnot(is.numeric(alpha))
	stopifnot(alpha>=0)
	if(nsim>0) testIC(nsim,alpha)

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
	
	if(cas==1) { #rectangle
		if(nsim==0) { #without CI
			res<-.C("corr_rect",
					as.integer(p$n),as.double(p$x),as.double(p$y),as.double(p$marks),
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
					as.integer(tmax),as.double(by),
					gm=double(tmax),km=double(tmax),
					PACKAGE="ads")
		}
		else { #with CI
			res<-.C("corr_rect_ic",
					as.integer(p$n),as.double(p$x),as.double(p$y),as.double(p$marks),
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
					as.integer(tmax),as.double(by),
					as.integer(nsim),as.double(alpha),
					gm=double(tmax),km=double(tmax),
					gmic1=double(tmax),gmic2=double(tmax),kmic1=double(tmax),kmic2=double(tmax),
					gmval=double(tmax),kmval=double(tmax),
					PACKAGE="ads")
		}
	}
	else if(cas==2) { #circle
		if(nsim==0) { #without CI
			res<-.C("corr_disq",
					as.integer(p$n),as.double(p$x),as.double(p$y),as.double(p$marks),
					as.double(x0),as.double(y0),as.double(r0),
					as.integer(tmax),as.double(by),
					gm=double(tmax),km=double(tmax),
					PACKAGE="ads")					
		}
		else { #with CI
			res<-.C("corr_disq_ic",
					as.integer(p$n),as.double(p$x),as.double(p$y),as.double(p$marks),
					as.double(x0),as.double(y0),as.double(r0),
					as.integer(tmax),as.double(by),
					as.integer(nsim),as.double(alpha),
					gm=double(tmax),km=double(tmax),
					gmic1=double(tmax),gmic2=double(tmax),kmic1=double(tmax),kmic2=double(tmax),
					gmval=double(tmax),kmval=double(tmax),
					PACKAGE="ads")
		}	
	}
	else if(cas==3) { #complex within rectangle
		if(nsim==0) { #without CI
			res<-.C("corr_tr_rect",
					as.integer(p$n),as.double(p$x),as.double(p$y),as.double(p$marks),
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by),
					gm=double(tmax),km=double(tmax),
					PACKAGE="ads")					
		}
		else { #with CI
			res<-.C("corr_tr_rect_ic",
					as.integer(p$n),as.double(p$x),as.double(p$y),as.double(p$marks),
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by),
					as.integer(nsim),as.double(alpha),
					gm=double(tmax),km=double(tmax),
					gmic1=double(tmax),gmic2=double(tmax),kmic1=double(tmax),kmic2=double(tmax),
					gmval=double(tmax),kmval=double(tmax),
					PACKAGE="ads")	
		}	
	}
	else if(cas==4) { #complex within circle
		if(nsim==0) { #without CI		
			res<-.C("corr_tr_disq",
					as.integer(p$n),as.double(p$x),as.double(p$y),as.double(p$marks),
					as.double(x0),as.double(y0),as.double(r0),
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by),
					gm=double(tmax),km=double(tmax),
					PACKAGE="ads")	
		}
		else { #with CI
			res<-.C("ripley_tr_disq_ic",
					as.integer(p$n),as.double(p$x),as.double(p$y),as.double(p$marks),
					as.double(x0),as.double(y0),as.double(r0),
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by),
					as.integer(nsim),as.double(alpha),
					gm=double(tmax),km=double(tmax),
					gmic1=double(tmax),gmic2=double(tmax),kmic1=double(tmax),kmic2=double(tmax),
					gmval=double(tmax),kmval=double(tmax),
					PACKAGE="ads")	
		}
	}
	# formatting results
	gm<-data.frame(obs=res$gm,theo=rep(0,tmax))
	km<-data.frame(obs=res$km,theo=rep(0,tmax))
	if(nsim>0) {
		gm<-cbind(gm,sup=res$gmic1,inf=res$gmic2,pval=res$gmval/(nsim+1))
		km<-cbind(km,sup=res$kmic1,inf=res$kmic2,pval=res$kmval/(nsim+1))
	}
	call<-match.call()
	res<-list(call=call,r=r,gm=gm,km=km)
	class(res)<-c("fads","kmfun")
	return(res)
}

ksfun<-function(p,upto,by,nsim=0,alpha=0.01) {
# checking for input parameters
	#options( CBoundsCheck = TRUE )
	stopifnot(inherits(p,"spp"))
	stopifnot(p$type=="multivariate")
	stopifnot(is.numeric(upto))
	stopifnot(upto>=1)
	stopifnot(is.numeric(by))
	stopifnot(by>0)
	r<-seq(by,upto,by)
	tmax<-length(r)
	stopifnot(is.numeric(nsim))
	stopifnot(nsim>=0)
	nsim<-testInteger(nsim)
	stopifnot(is.numeric(alpha))
	stopifnot(alpha>=0)
	if(nsim>0) testIC(nsim,alpha)
	
###faire test sur les marks
	
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
	intensity<-p$n/surface
	tabMarks<-levels(p$marks)
	nbMarks<-nlevels(p$marks)
#nbMarks<-length(tabMarks)
	marks<-as.numeric(p$marks)	
	freq<-as.vector(table(p$marks))
	D<-1-sum(freq*(freq-1))/(p$n*(p$n-1))

# computing Shimatani	
	if(cas==1) { #rectangle
		if(nsim==0) { #without CI
			res<-.C("shimatani_rect",
					as.integer(p$n),as.double(p$x),as.double(p$y),
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
					as.integer(tmax),as.double(by),
					as.integer(nbMarks),as.integer(marks),as.double(surface),gg=double(tmax),kk=double(tmax),erreur=integer(tmax),
					PACKAGE="ads")
		}
		else { #with CI
			res<-.C("shimatani_rect_ic",
					as.integer(p$n),as.double(p$x),as.double(p$y),as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
					as.integer(tmax),as.double(by), as.integer(nsim), as.double(alpha),
					as.integer(nbMarks),as.integer(marks),as.double(surface),as.double(D),
					gg=double(tmax),kk=double(tmax),gic1=double(tmax),gic2=double(tmax),kic1=double(tmax),kic2=double(tmax),
					gval=double(tmax),kval=double(tmax),erreur=integer(tmax),
					PACKAGE="ads")
		}
	}
	else if(cas==2) { #circle
		if(nsim==0) { #without CI
			res<-.C("shimatani_disq",
					as.integer(p$n),as.double(p$x),as.double(p$y),
					as.double(x0),as.double(y0),as.double(r0),
					as.integer(tmax),as.double(by),
					as.integer(nbMarks),as.integer(marks),as.double(surface),gg=double(tmax),kk=double(tmax),erreur=integer(tmax),
					PACKAGE="ads")
		}
		else { #with CI
			res<-.C("shimatani_disq_ic",
					as.integer(p$n),as.double(p$x),as.double(p$y),as.double(x0),as.double(y0),as.double(r0),
					as.integer(tmax),as.double(by), as.integer(nsim), as.double(alpha),
					as.integer(nbMarks),as.integer(marks),as.double(surface),as.double(D),
					gg=double(tmax),kk=double(tmax),gic1=double(tmax),gic2=double(tmax),kic1=double(tmax),kic2=double(tmax),
					gval=double(tmax),kval=double(tmax),erreur=integer(tmax),
					PACKAGE="ads")
		}
	}
	else if(cas==3) { #complex within rectangle
		if(nsim==0) { #without CI
			res<-.C("shimatani_tr_rect",
					as.integer(p$n),as.double(p$x),as.double(p$y),
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by),
					as.integer(nbMarks),as.integer(marks),as.double(surface),gg=double(tmax),kk=double(tmax),erreur=integer(tmax),
					PACKAGE="ads")
		}
		else { #with CI
			res<-.C("shimatani_tr_rect_ic",
					as.integer(p$n),as.double(p$x),as.double(p$y),as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by), as.integer(nsim), as.double(alpha),
					as.integer(nbMarks),as.integer(marks),as.double(surface),as.double(D),
					gg=double(tmax),kk=double(tmax),gic1=double(tmax),gic2=double(tmax),kic1=double(tmax),kic2=double(tmax),
					gval=double(tmax),kval=double(tmax),erreur=integer(tmax),
					PACKAGE="ads")
		}
	}
	else if(cas==4) { #complex within circle
		if(nsim==0) { #without CI		
			res<-.C("shimatani_tr_disq",
					as.integer(p$n),as.double(p$x),as.double(p$y),
					as.double(x0),as.double(y0),as.double(r0),
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by),
					as.integer(nbMarks),as.integer(marks),as.double(surface),gg=double(tmax),kk=double(tmax),erreur=integer(tmax),
					PACKAGE="ads")
		}
		else { #with CI
			res<-.C("shimatani_tr_disq_ic",
				   as.integer(p$n),as.double(p$x),as.double(p$y),as.double(x0),as.double(y0),as.double(r0),
				   as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
				   as.integer(tmax),as.double(by), as.integer(nsim), as.double(alpha),
				   as.integer(nbMarks),as.integer(marks),as.double(surface),as.double(D),
				   gg=double(tmax),kk=double(tmax),gic1=double(tmax),gic2=double(tmax),kic1=double(tmax),kic2=double(tmax),
				   gval=double(tmax),kval=double(tmax),erreur=integer(tmax),
				   PACKAGE="ads")
		}
	}		
	if(sum(res$erreur>0)){
		message("Error in ", appendLF=F)
		print(match.call())
		message("No neigbors within distance intervals:")
		print(paste(by*(res$erreur[res$erreur>0]-1),"-",by*res$erreur[res$erreur>0]))
		message("Increase argument 'by'")
		return(res=NULL)
	}
	gs<-data.frame(obs=res$gg/D,theo=rep(1,tmax))
	ks<-data.frame(obs=res$kk/D,theo=rep(1,tmax))
	if(nsim>0) {
		gs<-cbind(gs,sup=res$gic1/D,inf=res$gic2/D,pval=res$gval/(nsim+1))
		ks<-cbind(ks,sup=res$kic1/D,inf=res$kic2/D,pval=res$kval/(nsim+1))
	}
	call<-match.call()
	res<-list(call=call,r=r,gs=gs,ks=ks)
	class(res)<-c("fads","ksfun")
	return(res)
}


#################
#V2 that calls K12fun
##############
krfun<-function(p,upto,by,nsim=0,dis=NULL,H0=c("rl","se"),alpha=0.01) {
# checking for input parameters
	stopifnot(inherits(p,"spp"))
	stopifnot(p$type=="multivariate")
	stopifnot(is.numeric(upto))
	stopifnot(upto>=1)
	stopifnot(is.numeric(by))
	stopifnot(by>0)
	r<-seq(by,upto,by)
	tmax<-length(r)
	H0<-H0[1]
	stopifnot(H0=="se" || H0=="rl")
	ifelse(H0=="se",H0<-2,H0<-1)
	if(is.null(dis)) {
		stopifnot(H0==1)
		dis<-as.dist(matrix(1,nlevels(p$marks),nlevels(p$marks)))
		attr(dis,"Labels")<-levels(p$marks)
	}
	stopifnot(inherits(dis,"dist"))
	stopifnot(attr(dis,"Diag")==FALSE)
	stopifnot(attr(dis,"Upper")==FALSE)
	stopifnot(suppressWarnings(is.euclid(dis)))
###revoir tests sur dis	
	if(length(levels(p$marks))!=length(labels(dis))) {
		stopifnot(all(levels(p$marks)%in%labels(dis)))
#dis<-subsetdist(dis,which(labels(dis)%in%levels(p$marks)))
		dis<-subsetdist(dis,levels(p$marks))
		warning("matrix 'dis' have been subsetted to match with levels(p$marks)")
	}
#else if(any(labels(dis)!=levels(p$marks))) {
#		attr(dis,"Labels")<-levels(p$marks)
#		warning("levels(p$marks) have been assigned to attr(dis, ''Labels'')")
#	}
	stopifnot(is.numeric(nsim))
	stopifnot(nsim>=0)
	nsim<-testInteger(nsim)
	stopifnot(is.numeric(alpha))
	stopifnot(alpha>=0)
	if(nsim>0) testIC(nsim,alpha)
	
###faire test sur les marks
	
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
	intensity<-p$n/surface
	nbMarks<-nlevels(p$marks)
	marks<-as.numeric(p$marks) # => position du label dans levels(p$marks)
	dis<-as.dist(sortmat(dis,levels(p$marks)))
	HD<-suppressWarnings(divc(as.data.frame(unclass(table(p$marks))),sqrt(2*dis),scale=F)[1,1])
	HD<-HD*p$n/(p$n-1)
	dis<-as.vector(dis)
		
# computing Rao	
	if(cas==1) { #rectangle
		if(nsim==0) { #without CI
			res<-.C("rao_rect",
					as.integer(p$n),as.double(p$x),as.double(p$y),
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
					as.integer(tmax),as.double(by),as.integer(H0),
					as.integer(nbMarks),as.integer(marks),as.double(dis),as.double(surface),as.double(HD),
					gg=double(tmax),kk=double(tmax),gs=double(tmax),ks=double(tmax),erreur=integer(tmax),
					PACKAGE="ads")
		}
		else { #with CI
			res<-.C("rao_rect_ic",
					as.integer(p$n),as.double(p$x),as.double(p$y),as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
					as.integer(tmax),as.double(by), as.integer(nsim),as.integer(H0),as.double(alpha),
					as.integer(nbMarks),as.integer(marks),as.double(dis),as.double(surface),as.double(HD),
					gg=double(tmax),kk=double(tmax),gs=double(tmax),ks=double(tmax),gic1=double(tmax),gic2=double(tmax),kic1=double(tmax),kic2=double(tmax),
					gval=double(tmax),kval=double(tmax),erreur=integer(tmax),
					PACKAGE="ads")
		}
	}
	else if(cas==2) { #circle
		if(nsim==0) { #without CI
			res<-.C("rao_disq",
					as.integer(p$n),as.double(p$x),as.double(p$y),
					as.double(x0),as.double(y0),as.double(r0),
					as.integer(tmax),as.double(by),as.integer(H0),
					as.integer(nbMarks),as.integer(marks),as.double(dis),as.double(surface),as.double(HD),
					gg=double(tmax),kk=double(tmax),gs=double(tmax),ks=double(tmax),erreur=integer(tmax),
					PACKAGE="ads")
		}
		else { #with CI
			res<-.C("rao_disq_ic",
					as.integer(p$n),as.double(p$x),as.double(p$y),as.double(x0),as.double(y0),as.double(r0),
					as.integer(tmax),as.double(by), as.integer(nsim),as.integer(H0),as.double(alpha),
					as.integer(nbMarks),as.integer(marks),as.double(dis),as.double(surface),as.double(HD),
					gg=double(tmax),kk=double(tmax),gs=double(tmax),ks=double(tmax),gic1=double(tmax),gic2=double(tmax),kic1=double(tmax),kic2=double(tmax),
					gval=double(tmax),kval=double(tmax),erreur=integer(tmax),
					PACKAGE="ads")
		}
	}
	else if(cas==3) { #complex within rectangle
		if(nsim==0) { #without CI
			res<-.C("rao_tr_rect",
					as.integer(p$n),as.double(p$x),as.double(p$y),
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by),as.integer(H0),
					as.integer(nbMarks),as.integer(marks),as.double(dis),as.double(surface),as.double(HD),
					gg=double(tmax),kk=double(tmax),gs=double(tmax),ks=double(tmax),erreur=integer(tmax),
					PACKAGE="ads")
		}
		else { #with CI
			res<-.C("rao_tr_rect_ic",
					as.integer(p$n),as.double(p$x),as.double(p$y),
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by),as.integer(nsim), as.integer(H0),as.double(alpha),
					as.integer(nbMarks),as.integer(marks),as.double(dis),as.double(surface),as.double(HD),
					gg=double(tmax),kk=double(tmax),gs=double(tmax),ks=double(tmax),gic1=double(tmax),gic2=double(tmax),kic1=double(tmax),kic2=double(tmax),
					gval=double(tmax),kval=double(tmax),erreur=integer(tmax),
					PACKAGE="ads")
		}
	}
	else if(cas==4) { #complex within circle
		if(nsim==0) { #without CI	
			res<-.C("rao_tr_disq",
				   as.integer(p$n),as.double(p$x),as.double(p$y),
				   as.double(x0),as.double(y0),as.double(r0),
				   as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
				   as.integer(tmax),as.double(by),as.integer(H0),
				   as.integer(nbMarks),as.integer(marks),as.double(dis),as.double(surface),as.double(HD),
				   gg=double(tmax),kk=double(tmax),gs=double(tmax),ks=double(tmax),erreur=integer(tmax),
				   PACKAGE="ads")
		}
		else { #with CI (not based on K12)
			res<-.C("rao_tr_disq_ic",
					as.integer(p$n),as.double(p$x),as.double(p$y),
					as.double(x0),as.double(y0),as.double(r0),
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.integer(tmax),as.double(by),as.integer(nsim), as.integer(H0),as.double(alpha),
					as.integer(nbMarks),as.integer(marks),as.double(dis),as.double(surface),as.double(HD),
					gg=double(tmax),kk=double(tmax),gs=double(tmax),ks=double(tmax),gic1=double(tmax),gic2=double(tmax),kic1=double(tmax),kic2=double(tmax),
					gval=double(tmax),kval=double(tmax),serreur=integer(tmax),
					PACKAGE="ads")
		}
	}		
	if(sum(res$erreur>0)){
		message("Error in ", appendLF=F)
		print(match.call())
		message("No neigbors within distance intervals:")
		print(paste(by*(res$erreur[res$erreur>0]-1),"-",by*res$erreur[res$erreur>0]))
		message("Increase argument 'by'")
		return(res=NULL)
	}
	if(H0==1) {
		theog<-rep(1,tmax)
		theok<-rep(1,tmax)
	}
	if(H0==2) {
		theog<-res$gs
		theok<-res$ks
	}
	gr<-data.frame(obs=res$gg,theo=theog)
	kr<-data.frame(obs=res$kk,theo=theok)
	if(nsim>0) {
		gr<-cbind(gr,sup=res$gic1,inf=res$gic2,pval=res$gval/(nsim+1))
		kr<-cbind(kr,sup=res$kic1,inf=res$kic2,pval=res$kval/(nsim+1))
	}
	call<-match.call()
	res<-list(call=call,r=r,gr=gr,kr=kr)
	class(res)<-c("fads","krfun")
	return(res)
}

kdfun<-function(p,upto,by,dis,nsim=0,alpha=0.01) {
# checking for input parameters
	stopifnot(inherits(p,"spp"))
	stopifnot(p$type=="multivariate")
	if(min(p$x)<0)
		p$x<-p$x+abs(min(p$x))
	if(min(p$y)<0)
		p$y<-p$y+abs(min(p$y))
	stopifnot(is.numeric(upto))
	stopifnot(upto>=1)
	stopifnot(is.numeric(by))
	stopifnot(by>0)
	r<-seq(by,upto,by)
	tmax<-length(r)
	stopifnot(inherits(dis,"dist"))
	stopifnot(attr(dis,"Diag")==FALSE)
	stopifnot(attr(dis,"Upper")==FALSE)
	stopifnot(suppressWarnings(is.euclid(dis)))
###revoir tests sur dis	
	if(length(levels(p$marks))!=length(labels(dis))) {
		stopifnot(all(levels(p$marks)%in%labels(dis)))
#dis<-subsetdist(dis,which(labels(dis)%in%levels(p$marks)))
		dis<-subsetdist(dis,levels(p$marks))
		warning("matrix 'dis' have been subsetted to match with levels(p$marks)")
	}
#else if(any(labels(dis)!=levels(p$marks))) {
#		attr(dis,"Labels")<-levels(p$marks)
#		warning("levels(p$marks) have been assigned to attr(dis, ''Labels'')")
#	}
	stopifnot(is.numeric(nsim))
	stopifnot(nsim>=0)
	nsim<-testInteger(nsim)
	stopifnot(is.numeric(alpha))
	stopifnot(alpha>=0)
	if(nsim>0) testIC(nsim,alpha)
	
	surface<-area.swin(p$window)
	intensity<-p$n/surface
	nbMarks<-nlevels(p$marks)
	marks<-as.numeric(p$marks) # => position du label dans levels(p$marks)
	dis<-as.dist(sortmat(dis,levels(p$marks)))
	HD<-suppressWarnings(divc(as.data.frame(unclass(table(p$marks))),sqrt(2*dis),scale=F)[1,1])
	HD<-HD*p$n/(p$n-1)
	dis<-as.vector(dis)

###faire test sur les marks
	
	if(nsim==0) { #without CI
		res<-.C("shen",
			as.integer(p$n),as.double(p$x),as.double(p$y),
			as.integer(tmax),as.double(by),
			as.integer(nbMarks),as.integer(marks),as.double(dis),as.double(surface),as.double(HD),
			gd=double(tmax),kd=double(tmax),erreur=integer(tmax),
			PACKAGE="ads")
	}
	else { #with CI
		res<-.C("shen_ic",
			as.integer(p$n),as.double(p$x),as.double(p$y),
			as.integer(tmax),as.double(by), as.integer(nsim),as.double(alpha),
			as.integer(nbMarks),as.integer(marks),as.double(dis),as.double(surface),as.double(HD),
			gd=double(tmax),kd=double(tmax),gic1=double(tmax),gic2=double(tmax),kic1=double(tmax),kic2=double(tmax),
			gval=double(tmax),kval=double(tmax),erreur=integer(tmax),
			PACKAGE="ads")
	}	
	
	if(sum(res$erreur>0)){
		message("Error in ", appendLF=F)
		print(match.call())
		message("No neigbors within distance intervals:")
		print(paste(by*(res$erreur[res$erreur>0]-1),"-",by*res$erreur[res$erreur>0]))
		message("Increase argument 'by'")
		return(res=NULL)
	}
	gd<-data.frame(obs=res$gd,theo=rep(1,tmax))
	kd<-data.frame(obs=res$kd,theo=rep(1,tmax))
	if(nsim>0) {
		gd<-cbind(gd,sup=res$gic1,inf=res$gic2,pval=res$gval/(nsim+1))
		kd<-cbind(kd,sup=res$kic1,inf=res$kic2,pval=res$kval/(nsim+1))
	}
	call<-match.call()
	res<-list(call=call,r=r,gd=gd,kd=kd)
	class(res)<-c("fads","kdfun")
	return(res)
}
