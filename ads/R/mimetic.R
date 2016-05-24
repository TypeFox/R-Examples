#mimetic point process as in Goreaud et al. 2004
#RP 11/06/2013
###################################################

mimetic<-function(x,upto=NULL,by=NULL,prec=NULL,nsimax=3000,conv=50) {
# checking for input parameters 
	stopifnot(inherits(x,"fads")||inherits(x,"spp"))
	if(inherits(x,"fads")) {	
		call<-x$call
		p<-eval(call[[2]])
		upto<-call[[3]]
		by<-call[[4]]
		if(length(call)==6)
			prec<-call[[6]]
		else
			prec<-0.01
		lobs<-x$l$obs
		r<-x$r
		linit<-x
	}
	else if(inherits(x,"spp")) {
		p<-x
		if(is.null(prec))
			prec<-0.01
		else
			prec<-prec
		linit<-kfun(p=p,upto=upto,by=by,nsim=0,prec=prec)
		lobs<-linit$l$obs
		r<-linit$r
	}
	surface<-area.swin(p$window)
	tmax<-length(r)
#lobs<-lobs+r
	if("rectangle"%in%p$window$type) {
		xmin<-p$window$xmin
		xmax<-p$window$xmax
		ymin<-p$window$ymin
		ymax<-p$window$ymax
		stopifnot(upto<=(0.5*max((xmax-xmin),(ymax-ymin))))
		if ("complex"%in%p$window$type) {
			tri<-p$window$triangles
			nbTri<-nrow(tri)
			res<-.C("mimetic_tr_rect",
					as.integer(p$n),as.double(p$x),as.double(p$y),as.double(surface),
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.double(prec),as.integer(tmax),as.double(by),
					as.double(lobs),as.integer(nsimax),as.integer(conv),cost=double(nsimax),
					g=double(tmax),k=double(tmax),xx=double(p$n),yy=double(p$n),mess=as.integer(1),
					PACKAGE="ads")
		}
		else {
			res<-.C("mimetic_rect",
					as.integer(p$n),as.double(p$x),as.double(p$y),as.double(surface),
					as.double(xmin),as.double(xmax),as.double(ymin),as.double(ymax),
					as.double(prec),as.integer(tmax),as.double(by),
					as.double(lobs),as.integer(nsimax),as.integer(conv),cost=double(nsimax),
				g=double(tmax),k=double(tmax),xx=double(p$n),yy=double(p$n),mess=as.integer(1),
				PACKAGE="ads")
		}
	}
	else if("circle"%in%p$window$type) {
		x0<-p$window$x0
		y0<-p$window$y0
		r0<-p$window$r0
		stopifnot(upto<=r0)
		if ("complex"%in%p$window$type) {
			tri<-p$window$triangles
			nbTri<-nrow(tri)
			res<-.C("mimetic_tr_disq",
					as.integer(p$n),as.double(p$x),as.double(p$y),as.double(surface),
					as.double(x0),as.double(y0),as.double(r0),
					as.integer(nbTri),as.double(tri$ax),as.double(tri$ay),as.double(tri$bx),as.double(tri$by),as.double(tri$cx),as.double(tri$cy),
					as.double(prec),as.integer(tmax),as.double(by),
					as.double(lobs),as.integer(nsimax),as.integer(conv),cost=double(nsimax),
					g=double(tmax),k=double(tmax),xx=double(p$n),yy=double(p$n),mess=as.integer(1),
					PACKAGE="ads")
		}
		else {
			res<-.C("mimetic_disq",
					as.integer(p$n),as.double(p$x),as.double(p$y),as.double(surface),
					as.double(x0),as.double(y0),as.double(r0),as.double(prec),as.integer(tmax),as.double(by),
					as.double(lobs),as.integer(nsimax),as.integer(conv),cost=double(nsimax),
					g=double(tmax),k=double(tmax),xx=double(p$n),yy=double(p$n),mess=as.integer(1),
					PACKAGE="ads")	
		}
	}
	else
		stop("invalid window type")
	psim<-spp(res$x,res$y,window=p$window)
	lsim<-kfun(psim,upto,by,nsim=0,prec)
	cost<-res$cost
	call<-match.call()
	l<-data.frame(obs=linit$l$obs,sim=lsim$l$obs)
	fads<-list(r=lsim$r,l=l)
	class(fads)<-c("fads","mimetic")
	res<-list(call=call,fads=fads,spp=psim,cost=cost[cost>0])
	class(res)<-c("mimetic")
	return(res)
}

plot.mimetic<-function (x,cols,lty,main,sub,legend=TRUE,csize=1,cex.main=1.5,pos="top",...) {
	def.par <- par(no.readonly = TRUE)
	on.exit(par(def.par))
	mylayout<-layout(matrix(c(1,1,2,3,2,3),ncol=2,byrow=TRUE))
	main<-deparse(x$call,width.cutoff=100)
	plot.fads.mimetic(x$fads,main=main,cex.main=1.5*csize,pos=pos,...)
	plot(x$spp,main="x$spp (simulated)",cex.main=csize,...)
	barplot(x$cost,main=paste("x$cost (nsim=",length(x$cost),")",sep=""),cex.main=csize,...)
	
}