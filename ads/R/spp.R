spp<-function (x,y=NULL,window,triangles,marks,int2fac=TRUE) {
	if(is.list(x)) {
		stopifnot(length(x)==2)
		y<-x[[2]]
		x<-x[[1]]
	}
	else
		stopifnot(!is.null(y))
	stopifnot(is.numeric(x))
    stopifnot(is.numeric(y))
	stopifnot(length(x)==length(y))
	if(any(duplicated(cbind(x,y))))
		warning("duplicated (x,y) points")
#stopifnot(!duplicated(cbind(x,y)))
	if(!inherits(window,"swin")) {
		if(missing(triangles))
			w<-swin(window)
		else
			w<-swin(window,triangles)
	}
	else if("simple"%in%window$type&&!missing(triangles)) {
		if("rectangle"%in%window$type)
			w<-swin(c(window$xmin,window$ymin,window$xmax,window$ymax),triangles=triangles)
		else if("circle"%in%window$type)
			w<-swin(c(window$x0,window$y0,window$r0),triangles=triangles)
		else
			stop("invalid window type")
	}
	else
		w<-window
	stopifnot(any(ok<-inside.swin(x,y,w)!=FALSE))
	xout<-x[!ok]
	yout<-y[!ok]
	x<-x[ok]
	y<-y[ok]
	n<-length(x)
	nout<-length(xout)
	spp<-list(type="univariate",window=w,n=n,x=x,y=y)
	if(nout>0) spp<-c(spp,list(nout=nout,xout=xout,yout=yout))
	if(!missing(marks)) {
		stopifnot(length(marks)==(n+nout))
		if(!is.factor(marks)) {
			stopifnot(is.vector(marks))
			if(is.integer(marks)&&int2fac==TRUE) {
				marks<-as.factor(marks)
				warning("integer marks have been coerced to factor",call.=FALSE)
			}
			else if(is.character(marks))
				marks<-as.factor(marks)
			else {
				spp$type<-"marked"
				spp$marks <- marks[ok]
				if(nout>0)
					spp$marksout<-marks[!ok]
			}
		}
		if(is.factor(marks)) {
			spp$type<-"multivariate"
			names(marks)<-NULL
			spp$marks <- factor(marks[ok],exclude=NULL)
			if(nout>0)
				spp$marksout<-factor(marks[!ok],exclude=NULL)
		}
	}
	class(spp)<-"spp"
	return(spp)
}

print.spp<-function (x,...) {
	cat("Spatial point pattern:\n")
	str(x)
}

summary.spp<-function (object,...) {
    res<-list(type=object$type,window=summary(object$window))
	area<-res$window$area
	if(object$type=="multivariate")
		res$intensity<-summary(object$marks)/area
	else
		res$intensity<-object$n/area
	if(object$type=="marked")
		res$marks<-summary(object$marks)
	class(res)<-"summary.spp"
	return(res)
}

print.summary.spp<-function (x,...) {
	cat(paste("Spatial point pattern type:",x$type[1],"\n"))
	print(x$window)
	if(x$type=="multivariate") {
		cat("intensity:")
		print(array(x$intensity,dim=c(length(x$intensity),1),dimnames=list(paste(" .",names(x$intensity)),"")))
	}
	else
		cat(paste("intensity:",signif(x$intensity),"\n"))
	if(x$type=="marked") {
		cat("marks:\n")
		print(x$marks)
	}
}

plot.spp<-function (x,main,out=FALSE,use.marks=TRUE,cols,chars,cols.out,chars.out,maxsize,scale=TRUE,add=FALSE,legend=TRUE,csize=1,...) {
    stopifnot(x$n>0)
	if (missing(main)) 
        main <- deparse(substitute(x))
	#def.par <- par(no.readonly = TRUE)
	#on.exit(par(def.par))
	#if(options()$device=="windows")
	#	csize<-0.75*csize
	par(cex=csize)
	#print(par("cex"))
    if (!add) {
		if(out) {
			s<-summary(x$window)
			e<-max(c(diff(c(min(x$xout),s$xrange[1])),diff(c(s$xrange[2],max(x$xout))),
			diff(c(min(x$yout),s$yrange[1])),diff(c(s$yrange[2],max(x$yout)))))
			plot.swin(x$window,main,e,scale,csize=csize,...)
		}
		else if(x$type!="univariate"&&legend==TRUE) {
			s<-summary(x$window)
			e<-0.2*(s$yrange[2]-s$yrange[1])
			plot.swin(x$window,main,e,scale,csize=csize,...)
		}
		else
			plot.swin(x$window,main,scale=scale,csize=csize,...)
	}
	if(x$type=="univariate"||!use.marks) {
		if(!missing(cols))
			stopifnot(length(cols)==1)
		if(!missing(chars))
			stopifnot(length(chars)==1)
		points(x$x,x$y,pch=ifelse(missing(chars),1,chars),col=ifelse(missing(cols),"black",cols),cex=ifelse(missing(maxsize),1,maxsize),...)
		if(out) {
			if(!missing(cols.out))
				stopifnot(length(cols.out)==1)
			if(!missing(chars.out))
				stopifnot(length(chars.out)==1)
			points(x$xout,x$yout,pch=ifelse(missing(chars.out),2,chars.out),col=ifelse(missing(cols.out),"red",cols.out),
				cex=ifelse(missing(maxsize),1,maxsize),...)
		}
    }
	else if(x$type=="multivariate") {
		if(!missing(cols))
			stopifnot(length(cols)==nlevels(x$marks))
		if(!missing(chars))
			stopifnot(length(chars)==nlevels(x$marks))
		for(i in 1:nlevels(x$marks)) {
			rel<-(x$marks==levels(x$marks)[i])
			points(x$x[rel],x$y[rel],pch=ifelse(missing(chars),i,chars[i]),col=ifelse(missing(cols),i,cols[i]),
				cex=ifelse(missing(maxsize),1,maxsize),...)
		}
		if(out) {
			if(!missing(cols.out))
				stopifnot(length(cols.out)==nlevels(x$marksout))
			if(!missing(chars.out))
				stopifnot(length(chars.out)==nlevels(x$marksout))
			for(i in 1:nlevels(x$marksout)) {
				rel<-(x$marksout==levels(x$marksout)[i])
				points(x$xout[rel],x$yout[rel],pch=ifelse(missing(chars.out),i,chars.out[i]),
					col=ifelse(missing(cols.out),i,cols.out[i]),cex=ifelse(missing(maxsize),1,maxsize),...)
			}
		}
		if(legend) {
			xl<-c(s$xrange[1],s$xrange[2])
			yl<-c(s$yrange[2]*1.15,s$yrange[2])
			#xl<-c(p$window$xmin,p$window$xmax)
			#yl<-c(p$window$ymax*1.1,p$window$ymax)
			legend(x=xl,y=yl,levels(x$marks),cex=1,bty="n",pch=if(missing(chars)) c(1:nlevels(x$marks)) 
				else chars,col=if(missing(cols)) c(1:nlevels(x$marks)) else cols,horiz=TRUE,xjust=0.5,...)
		}
	}
	else if(x$type=="marked") {
		if(!missing(cols))
			stopifnot(length(cols)==1)
		if(!missing(chars))
			stopifnot(length(chars)==1)
		if(out) {
			#ms<-adjust.marks.size(c(x$marks,x$marksout),x$window,if(!missing(maxsize)) maxsize)
			ms<-adjust.marks.size(c(x$marks,x$marksout),x$window)
			if(!missing(maxsize))
				ms<-ms*maxsize
			msout<-ms[(x$n+1):(x$n+x$nout)]
			ms<-ms[1:x$n]
		}	
		else {
			#ms<-adjust.marks.size(x$marks,x$window,if(!missing(maxsize)) maxsize)
			ms<-adjust.marks.size(x$marks,x$window)
			if(!missing(maxsize))
				ms<-ms*maxsize
		}
		neg<-(x$marks<0)
		if(missing(chars)||chars=="circles") {
			if(any(neg))
				symbols(x$x[neg],x$y[neg],circles=-ms[neg]/2,fg=ifelse(missing(cols),1,cols),inches=FALSE,add=TRUE,...)
			if(any(!neg))
				symbols(x$x[!neg],x$y[!neg],circles=ms[!neg]/2,fg=ifelse(missing(cols),1,cols),bg=ifelse(missing(cols),1,cols),inches=FALSE,add=TRUE,...)
		}
		else if(chars=="squares") {
			if(any(neg))
				symbols(x$x[neg],x$y[neg],squares=-ms[neg],fg=ifelse(missing(cols),1,cols),inches=FALSE,add=TRUE,...)
			if(any(!neg))
				symbols(x$x[!neg],x$y[!neg],squares=ms[!neg],fg=ifelse(missing(cols),1,cols),bg=ifelse(missing(cols),1,cols),inches=FALSE,add=TRUE,...)
		}
		else
			stopifnot(chars%in%c("circles","squares"))
		if(out) {
			neg<-(x$marksout<0)
			if(missing(chars.out)||chars.out=="circles") {
				if(any(neg))
					symbols(x$xout[neg],x$yout[neg],circles=-msout[neg]/2,fg=ifelse(missing(cols.out),2,cols.out),inches=FALSE,add=TRUE,...)
				if(any(!neg))
					symbols(x$xout[!neg],x$yout[!neg],circles=msout[!neg]/2,fg=ifelse(missing(cols.out),2,cols.out),
					bg=ifelse(missing(cols.out),2,cols.out),inches=FALSE,add=TRUE,...)

			}
			else if(chars.out=="squares") {
				if(any(neg))
					symbols(x$xout[neg],x$yout[neg],squares=-msout[neg],fg=ifelse(missing(cols.out),2,cols.out),inches=FALSE,add=TRUE,...)
				if(any(!neg))
					symbols(x$xout[!neg],x$yout[!neg],squares=msout[!neg],fg=ifelse(missing(cols.out),2,cols.out),
					bg=ifelse(missing(cols.out),2,cols.out),inches=FALSE,add=TRUE,...)

			}
			else
				stopifnot(chars.out%in%c("circles","squares"))
		}
		if(legend) {
			mid<-(s$xrange[2]-s$xrange[1])/2
			xl<-c(mid-0.5*mid,mid,mid+0.5*mid)
			yl<-rep((s$yrange[2]*1.15),3)
			lm<-range(abs(x$marks))
			lm<-c(lm[1],mean(lm),lm[2])
			lms<-range(abs(ms))
			lms<-c(lms[1],mean(lms),lms[2])
			if(missing(chars)||chars=="circles") {
				symbols(xl,yl,circles=lms/2,fg=ifelse(missing(cols),1,cols),bg=ifelse(missing(cols),1,cols),inches=FALSE,add=TRUE,...)
				text(c(xl[1]+lms[1]+1,xl[2]+lms[2]+1,xl[3]+lms[3]+1),yl,labels=signif(lm,2),pos=4,cex=1)
				if(any(neg)) {
					symbols(xl,yl*0.93,circles=lms/2,fg=ifelse(missing(cols),1,cols),inches=FALSE,add=TRUE,...)
					text(c(xl[1]+lms[1],xl[2]+lms[2],xl[3]+lms[3]),yl*0.93,labels=signif(-lm,2),pos=4,cex=1)
				}
			}
			else if(chars=="squares") {
				symbols(xl,yl,squares=lms,fg=ifelse(missing(cols),1,cols),bg=ifelse(missing(cols),1,cols),inches=FALSE,add=TRUE,...)
				text(c(xl[1]+lms[1]+1,xl[2]+lms[2]+1,xl[3]+lms[3]+1),yl,labels=signif(lm,2),pos=4,cex=1)
				if(any(neg)) {
					symbols(xl,yl*0.93,squares=lms,fg=ifelse(missing(cols),1,cols),inches=FALSE,add=TRUE,...)
					text(c(xl[1]+lms[1],xl[2]+lms[2],xl[3]+lms[3]),yl*0.93,labels=signif(-lm,2),pos=4,cex=1)
				}
			}
		}

	}
	else
		stop("invalid point pattern type")
}

ppp2spp<-function(p) {
	stopifnot(inherits(p,"ppp"))
	sw<-owin2swin(p$window)
	if(is.null(p$marks))
		sp<-spp(p$x,p$y,sw)
	else
		sp<-spp(p$x,p$y,window=sw,marks=p$marks)
	return(sp)
}
