plot.vads<-function(x,main,opt,select,chars,cols,maxsize,char0,col0,legend,csize,...) {
	UseMethod("plot.vads")
}

plot.vads.dval<-function (x,main,opt=c("dval","cval"),select,chars=c("circles","squares"),cols,maxsize,char0,col0,legend=TRUE,csize=1,...) {
	if(!missing(select)) {
		d<-c()
		for(i in 1:length(select)) {
			select.in.r<-c()
			for(j in 1:length(x$r)) {
				select.in.r<-c(select.in.r,ti<-isTRUE(all.equal(select[i],x$r[j])))
				if(ti)
					d<-c(d,j)
			}
			stopifnot(any(select.in.r==TRUE))
		}
	}	
	else
		d<-rank(x$r)
	nd<-length(d)
	nf<-ceiling(sqrt(nd))
	stopifnot(opt%in%c("dval","cval"))
	opt<-opt[1]
	stopifnot(chars%in%c("circles","squares"))
	chars<-chars[1]
	ifelse(opt=="dval",val<-x$dval[,d],val<-x$cval[,d])
	v<-val
	val<-data.frame(adjust.marks.size(val,x$window,if(!missing(maxsize)) maxsize))
	def.par <- par(no.readonly = TRUE)
	on.exit(par(def.par))
	#if(options()$device=="windows")
	#	csize<-0.75*csize
	if (missing(main)) 
        main <- deparse(substitute(x))
	mylayout<-layout(matrix(c(rep(1,nf),seq(2,((nf*nf)+1),1)),(nf+1),nf,byrow=TRUE))
	s<-summary(x$window)
	par(mar=c(0.1,0.1,1,0.1),cex=csize)
	plot(s$xrange,s$yrange,type="n",axes=FALSE,asp=1/nf)
	legend("center","",cex=1.5,bty="n",horiz=TRUE,title=main,...)
	if(legend) {
			mid<-(s$xrange[2]-s$xrange[1])/2
			xl<-c(mid-0.5*mid,mid,mid+0.5*mid)
			yl<-rep(s$xrange[2]*0.25,3)
			lm<-range(v[v>0])
			lm<-c(lm[1],mean(lm),lm[2])
			lms<-range(val[val>0])
			lms<-c(lms[1],mean(lms),lms[2])
			if(missing(chars)||chars=="circles") {
				symbols(xl,yl,circles=sqrt(lms),fg=ifelse(missing(cols),1,cols),bg=ifelse(missing(cols),1,cols),inches=FALSE,add=TRUE,...)
				text(c(xl[1]+lms[1],xl[2]+lms[2],xl[3]+lms[3]),yl,labels=signif(lm,2),pos=4,cex=1.5)
			}
			else if(chars=="squares") {
				symbols(xl,yl,squares=1.5*sqrt(lms),fg=ifelse(missing(cols),1,cols),bg=ifelse(missing(cols),1,cols),inches=FALSE,add=TRUE,...)
				text(c(xl[1]+lms[1],xl[2]+lms[2],xl[3]+lms[3]),yl,labels=signif(lm,2),pos=4,cex=1.5)
			}
	}
	#if(!is.null(main)) {
	#	mylayout<-layout(matrix(c(rep(1,nf),seq(2,((nf*nf)+1),1)),(nf+1),nf,byrow=TRUE))
	#	plot.default(x$xy$x,x$xy$y,type="n",axes=FALSE)
	#	text(mean(range(x$xy$x)),mean(range(x$xy$y)),pos=3,cex=2,labels=main)
		#ajouter une lŽgende ???
	#}
	#else
	#	mylayout<-layout(matrix(seq(1,(nf*nf),1),nf,nf,byrow=TRUE))
	ifelse(missing(cols),cols<-1,cols<-cols[1])
	if(!missing(char0)||!missing(col0)) {
		ifelse(missing(col0),col0<-cols,col0<-col0[1])	
		if(missing(char0))
			char0<-3
	}
	for(i in 1:nd) {
		plot.swin(x$window,main=paste("r =",x$r[d[i]]),scale=FALSE,csize=0.66*csize,...)
		nort<-(val[,i]==0)
		if(!missing(char0)&&any(nort))
			points(x$xy$x[nort],x$xy$y[nort],pch=char0,col=col0,...)
		if(any(!nort)) {
			if(chars=="circles")
					symbols(x$xy$x[!nort],x$xy$y[!nort],circles=nf*sqrt(val[!nort,i]),
					fg=cols,bg=cols,inches=FALSE,add=TRUE,...)
			else if(chars=="squares")
					symbols(x$xy$x[!nort],x$xy$y[!nort],squares=1.5*nf*sqrt(val[!nort,i]),
					fg=cols,bg=cols,inches=FALSE,add=TRUE,...)
		}
	}
	## mŽthode en courbes de niveaux ?
}

plot.vads.kval<-function (x,main,opt=c("lval","kval","nval","gval"),select,chars=c("circles","squares"),cols,maxsize,char0,col0,legend=TRUE,csize=1,...) {
	if(!missing(select)) {
		d<-c()
		for(i in 1:length(select)) {
			select.in.r<-c()
			for(j in 1:length(x$r)) {
				select.in.r<-c(select.in.r,ti<-isTRUE(all.equal(select[i],x$r[j])))
				if(ti)
					d<-c(d,j)
			}
			stopifnot(any(select.in.r==TRUE))
		}
	}	
	else
		d<-rank(x$r)
	nd<-length(d)
	nf<-ceiling(sqrt(nd))
	opt<-opt[1]
	stopifnot(chars%in%c("circles","squares"))
	chars<-chars[1]
	
	if(opt=="lval")
		val<-x$lval[,d]
	else if(opt=="kval")
		val<-x$kval[,d]
	else if(opt=="nval")
		val<-x$nval[,d]
	else if(opt=="gval")
		val<-x$gval[,d]
	else
		stopifnot(opt%in%c("lval","kval","nval","gval"))
	v<-val
	#val<-data.frame(adjust.marks.size(val,x$window,if(!missing(maxsize)) maxsize))
	val<-data.frame(adjust.marks.size(val,x$window))
	if(!missing(maxsize))
		val<-val*maxsize
	def.par <- par(no.readonly = TRUE)
	on.exit(par(def.par))
	#if(options()$device=="windows")
	#	csize<-0.75*csize
	if (missing(main)) 
        main <- deparse(substitute(x))
	mylayout<-layout(matrix(c(rep(1,nf),seq(2,((nf*nf)+1),1)),(nf+1),nf,byrow=TRUE))
	s<-summary(x$window)
	par(mar=c(0.1,0.1,1,0.1),cex=csize)
	plot.default(s$xrange,s$yrange,type="n",axes=FALSE,asp=1/nf)
	legend("center","",cex=1.5,bty="n",horiz=TRUE,title=main,...)
	if(legend) {
			mid<-(s$xrange[2]-s$xrange[1])/2
			xl<-c(mid-0.5*mid,mid,mid+0.5*mid)
			yl<-rep(s$xrange[2]*0.25,3)
			lm<-range(abs(v)[abs(v)>0])
			lm<-c(lm[1],mean(lm),lm[2])
			lms<-range(abs(val)[abs(val)>0])
			lms<-c(lms[1],mean(lms),lms[2])
			if(missing(chars)||chars=="circles") {
				symbols(xl,yl,circles=sqrt(lms),fg=ifelse(missing(cols),1,cols),bg=ifelse(missing(cols),1,cols),inches=FALSE,add=TRUE,...)
				text(c(xl[1]+lms[1]+1,xl[2]+lms[2]+1,xl[3]+lms[3]+1),yl,labels=signif(lm,2),pos=4,cex=1)
				symbols(xl,yl*0.5,circles=sqrt(lms),fg=ifelse(missing(cols),1,cols),bg="white",inches=FALSE,add=TRUE,...)
				text(c(xl[1]+lms[1],xl[2]+lms[2],xl[3]+lms[3]),yl*0.5,labels=signif(-lm,2),pos=4,cex=1)

			}
			else if(chars=="squares") {
				symbols(xl,yl,squares=1.5*sqrt(lms),fg=ifelse(missing(cols),1,cols),bg=ifelse(missing(cols),1,cols),inches=FALSE,add=TRUE,...)
				text(c(xl[1]+lms[1]+1,xl[2]+lms[2]+1,xl[3]+lms[3]+1),yl,labels=signif(lm,2),pos=4,cex=1)
				symbols(xl,yl*0.5,squares=1.5*sqrt(lms),fg=ifelse(missing(cols),1,cols),bg="white",inches=FALSE,add=TRUE,...)
				text(c(xl[1]+lms[1],xl[2]+lms[2],xl[3]+lms[3]),yl*0.5,labels=signif(-lm,2),pos=4,cex=1)

			}
	}
	#if(!is.null(main)) {
	#	mylayout<-layout(matrix(c(rep(1,nf),seq(2,((nf*nf)+1),1)),(nf+1),nf,byrow=TRUE))
	#	plot.default(x$xy$x,x$xy$y,type="n",axes=FALSE)
	#	text(mean(range(x$xy$x)),mean(range(x$xy$y)),pos=3,cex=2,labels=main)
		#ajouter une lŽgende ???
	#}
	#else
	#	mylayout<-layout(matrix(seq(1,(nf*nf),1),nf,nf,byrow=TRUE))
	ifelse(missing(cols),cols<-1,cols<-cols[1])
	if(!missing(char0)||!missing(col0)) {
		ifelse(missing(col0),col0<-cols,col0<-col0[1])	
		if(missing(char0))
			char0<-3
	}				
	for(i in 1:nd) {
		plot.swin(x$window,main=paste("r =",x$r[d[i]]),scale=FALSE,csize=0.66*csize,...)
		nort<-(val[,i]==0)
		neg<-(val[,i]<0)
		if(!missing(char0)&&any(nort))
				points(x$xy$x[nort],x$xy$y[nort],pch=char0,col=col0,...)
		if(any(!nort)) {		
			if(chars=="circles") {
				if(any(!neg))
					symbols(x$xy$x[(!neg&!nort)],x$xy$y[(!neg&!nort)],circles=nf*sqrt(abs(val[(!neg&!nort),i])),
					fg=cols,bg=cols,inches=FALSE,add=TRUE,...)
				if(any(neg))
					symbols(x$xy$x[(neg&!nort)],x$xy$y[(neg&!nort)],circles=nf*sqrt(abs(val[(neg&!nort),i])),
					fg=cols,bg="white",inches=FALSE,add=TRUE,...)
			}
			else if(chars=="squares") {
				if(any(!neg))
					symbols(x$xy$x[(!neg&!nort)],x$xy$y[(!neg&!nort)],squares=1.5*nf*sqrt(abs(val[(!neg&!nort),i])),
					fg=cols,bg=cols,inches=FALSE,add=TRUE,...)
				if(any(neg))
					symbols(x$xy$x[(neg&!nort)],x$xy$y[(neg&!nort)],squares=1.5*nf*sqrt(abs(val[(neg&!nort),i])),
					fg=cols,bg="white",inches=FALSE,add=TRUE,...)
			}	
		}
	}
	## mŽthode en courbes de niveaux ?
}

plot.vads.k12val<-function (x,main,opt=c("lval","kval","nval","gval"),select,chars=c("circles","squares"),cols,maxsize,char0,col0,legend=TRUE,csize=1,...) {
	if(!missing(select)) {
		d<-c()
		for(i in 1:length(select)) {
			select.in.r<-c()
			for(j in 1:length(x$r)) {
				select.in.r<-c(select.in.r,ti<-isTRUE(all.equal(select[i],x$r[j])))
				if(ti)
					d<-c(d,j)
			}
			stopifnot(any(select.in.r==TRUE))
		}
	}	
	else
		d<-rank(x$r)
	nd<-length(d)
	nf<-ceiling(sqrt(nd))
	opt<-opt[1]
	stopifnot(chars%in%c("circles","squares"))
	chars<-chars[1]
	
	if(opt=="lval")
		val<-x$l12val[,d]
	else if(opt=="kval")
		val<-x$k12val[,d]
	else if(opt=="nval")
		val<-x$n12val[,d]
	else if(opt=="gval")
		val<-x$g12val[,d]
	else
		stopifnot(opt%in%c("lval","kval","nval","gval"))
	v<-val
	#val<-data.frame(adjust.marks.size(val,x$window,if(!missing(maxsize)) maxsize))
	val<-data.frame(adjust.marks.size(val,x$window))
	if(!missing(maxsize))
		val<-val*maxsize
	def.par <- par(no.readonly = TRUE)
	on.exit(par(def.par))
	#if(options()$device=="windows")
	#	csize<-0.75*csize
	if (missing(main)) 
        main <- deparse(substitute(x))
	mylayout<-layout(matrix(c(rep(1,nf),seq(2,((nf*nf)+1),1)),(nf+1),nf,byrow=TRUE))
	s<-summary(x$window)
	par(mar=c(0.1,0.1,1,0.1),cex=csize)
	plot.default(s$xrange,s$yrange,type="n",axes=FALSE,asp=1/nf)
	legend("center","",cex=1.5,bty="n",horiz=TRUE,title=main,...)
	if(legend) {
			mid<-(s$xrange[2]-s$xrange[1])/2
			xl<-c(mid-0.5*mid,mid,mid+0.5*mid)
			yl<-rep(s$xrange[2]*0.25,3)
			lm<-range(abs(v)[abs(v)>0])
			lm<-c(lm[1],mean(lm),lm[2])
			lms<-range(abs(val)[abs(val)>0])
			lms<-c(lms[1],mean(lms),lms[2])
			if(missing(chars)||chars=="circles") {
				symbols(xl,yl,circles=sqrt(lms),fg=ifelse(missing(cols),1,cols),bg=ifelse(missing(cols),1,cols),inches=FALSE,add=TRUE,...)
				text(c(xl[1]+lms[1]+1,xl[2]+lms[2]+1,xl[3]+lms[3]+1),yl,labels=signif(lm,2),pos=4,cex=1)
				symbols(xl,yl*0.5,circles=sqrt(lms),fg=ifelse(missing(cols),1,cols),bg="white",inches=FALSE,add=TRUE,...)
				text(c(xl[1]+lms[1],xl[2]+lms[2],xl[3]+lms[3]),yl*0.5,labels=signif(-lm,2),pos=4,cex=1)
			}
			else if(chars=="squares") {
				symbols(xl,yl,squares=1.5*sqrt(lms),fg=ifelse(missing(cols),1,cols),bg=ifelse(missing(cols),1,cols),inches=FALSE,add=TRUE,...)
				text(c(xl[1]+lms[1],xl[2]+lms[2],xl[3]+lms[3]),yl,labels=signif(lm,2),pos=4,cex=1)
				symbols(xl,yl*0.5,squares=1.5*sqrt(lms),fg=ifelse(missing(cols),1,cols),bg="white",inches=FALSE,add=TRUE,...)
				text(c(xl[1]+lms[1],xl[2]+lms[2],xl[3]+lms[3]),yl*0.5,labels=signif(-lm,2),pos=4,cex=1)
			}
	}
	#if(!is.null(main)) {
	#	mylayout<-layout(matrix(c(rep(1,nf),seq(2,((nf*nf)+1),1)),(nf+1),nf,byrow=TRUE))
	#	plot.default(x$xy$x,x$xy$y,type="n",axes=FALSE)
	#	text(mean(range(x$xy$x)),mean(range(x$xy$y)),pos=3,cex=2,labels=main)
		#ajouter une lŽgende ???
	#}
	#else
	#	mylayout<-layout(matrix(seq(1,(nf*nf),1),nf,nf,byrow=TRUE))
	ifelse(missing(cols),cols<-1,cols<-cols[1])
	if(!missing(char0)||!missing(col0)) {
		ifelse(missing(col0),col0<-cols,col0<-col0[1])	
		if(missing(char0))
			char0<-3
	}			
	for(i in 1:nd) {
		plot.swin(x$window,main=paste("r =",x$r[d[i]]),scale=FALSE,csize=0.66*csize,...)
		nort<-(val[,i]==0)
		neg<-(val[,i]<0)
		if(!missing(char0)&&any(nort))
				points(x$xy$x[nort],x$xy$y[nort],pch=char0,col=col0,...)
		if(any(!nort)) {		
			if(chars=="circles") {
				if(any(!neg))
					symbols(x$xy$x[(!neg&!nort)],x$xy$y[(!neg&!nort)],circles=nf*sqrt(abs(val[(!neg&!nort),i])),
					fg=cols,bg=cols,inches=FALSE,add=TRUE,...)
				if(any(neg))
					symbols(x$xy$x[(neg&!nort)],x$xy$y[(neg&!nort)],circles=nf*sqrt(abs(val[(neg&!nort),i])),
					fg=cols,bg="white",inches=FALSE,add=TRUE,...)
			}
			else if(chars=="squares") {
				if(any(!neg))
					symbols(x$xy$x[(!neg&!nort)],x$xy$y[(!neg&!nort)],squares=1.5*nf*sqrt(abs(val[(!neg&!nort),i])),
					fg=cols,bg=cols,inches=FALSE,add=TRUE,...)
				if(any(neg))
					symbols(x$xy$x[(neg&!nort)],x$xy$y[(neg&!nort)],squares=1.5*nf*sqrt(abs(val[(neg&!nort),i])),
					fg=cols,bg="white",inches=FALSE,add=TRUE,...)
			}	
		}
	}
	## mŽthode en courbes de niveaux ?
}
