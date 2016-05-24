plot.fads<-function (x,opt,cols,lty,main,sub,legend,csize,...) {
	UseMethod("plot.fads")
}

plot.fads.kfun<-function (x,opt=c("all","L","K","n","g"),cols,lty,main,sub,legend=TRUE,csize=1,...) {
	ifelse(!is.null(x$call$nsim)&&(x$call$nsim>0),ci<-TRUE,ci<-FALSE)
	def.par <- par(no.readonly = TRUE)
	on.exit(par(def.par))
	#if(options()$device=="windows")
	#	csize<-0.75*csize
	opt<-opt[1]
	if(opt=="all")
		mylayout<-layout(matrix(c(1,1,1,1,2,2,3,3,2,2,3,3,4,4,5,5,4,4,5,5),ncol=4,byrow=TRUE))
	else if(opt%in%c("L","K","n","g"))
		mylayout<-layout(matrix(c(1,1,1,1,rep(2,16)),ncol=4,byrow=TRUE))
	else
		stopifnot(opt%in%c("all","L","K","n","g"))
	if(missing(cols))
		cols=c(1,2,3)
	else if(length(cols)!=3)
		cols=c(cols,cols,cols)
	if(missing(lty))
			lty=c(1,3,2)
	else if(length(lty)!=3)
		lty=c(lty,lty,lty)
	if(missing(main))
		main<-deparse(x$call,width.cutoff=100)
	if(missing(sub))
		sub<-c("pair density function","second-order neighbour density function","Ripley's K-function","L-function : sqrt[K(r)/pi]-r")
	if(ci) {
		alpha<-x$call[["alpha"]]
		p<-ifelse(!is.null(alpha),signif(100*(1-alpha),digits=6),99)
		par(mar=c(0.1,0.1,0.1,0.1),cex=csize)
		plot(x$r,x$g$obs/2,type="n",axes=FALSE,xlab="",ylab="")
		if(legend)
			legend("center",c("obs","theo (CSR)",paste(p,"% CI of CSR")),cex=1.5,lty=lty[1:3],bty="n",horiz=TRUE,title=main,col=cols[1:3],...)
		else
			legend("center","",cex=1.5,bty="n",horiz=TRUE,title=main,...)
		par(mar=c(5,5,0.1,2),cex=ifelse(opt%in%c("all"),0.75*csize,csize))
		if(opt%in%c("all","g")) { # g-function
			lim<-range(x$g[,1:4])
			plot(x$r,x$g$obs,ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[1]),type="n",xlab="distance (r)",ylab="g(r)",cex.lab=1.25,...)
			lines(x$r,x$g$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$g$theo,lty=lty[2],col=cols[2],...)
			lines(x$r,x$g$sup,lty=lty[3],col=cols[3],...)
			lines(x$r,x$g$inf,lty=lty[3],col=cols[3],...)	
		}
		if(opt%in%c("all","n")) {# n-function
			lim<-range(x$n[,1:4])
			plot(x$r,x$n$obs,ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[2]),type="n",xlab="distance (r)",ylab="n(r)",cex.lab=1.25,...)
			lines(x$r,x$n$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$n$theo,lty=lty[2],col=cols[2],...)
			lines(x$r,x$n$sup,lty=lty[3],col=cols[3],...)
			lines(x$r,x$n$inf,lty=lty[3],col=cols[3],...)
		}
		if(opt%in%c("all","K")) { # K-function
			plot(x$r,x$k$obs,ylim=range(x$k[,1:4]),main=paste("\n\n",sub[3]),type="n",xlab="distance (r)",ylab="K(r)",cex.lab=1.25,...)
			lines(x$r,x$k$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$k$theo,lty=lty[2],col=cols[2],...)
			lines(x$r,x$k$sup,lty=lty[3],col=cols[3],...)
			lines(x$r,x$k$inf,lty=lty[3],col=cols[3],...)
		}
		if(opt%in%c("all","L")) { # L-function
			lim<-range(x$l[,1:4])
			plot(x$r,x$l$obs,ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[4]),type="n",xlab="distance (r)",ylab="L(r)",cex.lab=1.25,...)
			lines(x$r,x$l$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$l$theo,lty=lty[2],col=cols[2],...)
			lines(x$r,x$l$sup,lty=lty[3],col=cols[3],...)
			lines(x$r,x$l$inf,lty=lty[3],col=cols[3],...)
		}
	}
	else {
		par(mar=c(0.1,0.1,0.1,0.1),cex=csize)
		plot(x$r,x$g$obs/2,type="n",axes=FALSE,xlab="",ylab="")
		if(legend)
			legend("center",c("obs","theo (CSR)"),cex=1.5,lty=lty[1:2],bty="n",horiz=TRUE,title=main,col=cols[1:2],...)
		else
			legend("center","",cex=1.5,bty="n",horiz=TRUE,title=main,...)
		par(mar=c(5,5,0.1,2),cex=ifelse(opt%in%c("all"),0.75*csize,csize))
		if(opt%in%c("all","g")) { # g-function
			lim<-range(x$g)
			plot(x$r,x$g$obs,ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[1]),type="n",xlab="distance (r)",ylab="g(r)",cex.lab=1.25,...)
			lines(x$r,x$g$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$g$theo,lty=lty[2],col=cols[2],...)
		}
		if(opt%in%c("all","n")) { # n-function
			lim<-range(x$n)
			plot(x$r,x$n$obs,ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[2]),type="n",xlab="distance (r)",ylab="n(r)",cex.lab=1.25,...)
			lines(x$r,x$n$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$n$theo,lty=lty[2],col=cols[2],...)
		}
		if(opt%in%c("all","K")) { # k-function
			plot(x$r,x$k$obs,ylim=range(x$k),main=paste("\n\n",sub[3]),type="n",xlab="distance (r)",ylab="K(r)",cex.lab=1.25,...)
			lines(x$r,x$k$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$k$theo,lty=lty[2],col=cols[2],...)
		}
		if(opt%in%c("all","L")) { # L-function
			lim<-range(x$l)
			plot(x$r,x$l$obs,ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[4]),type="n",xlab="distance (r)",ylab="L(r)",cex.lab=1.25,...)
			lines(x$r,x$l$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$l$theo,lty=lty[2],col=cols[2],...)
		}
	}	
}

plot.fads.k12fun<-function(x,opt=c("all","L","K","n","g"),cols,lty,main,sub,legend=TRUE,csize=1,...) {
	#ifelse(!is.null(x$call[["nsim"]]),ci<-TRUE,ci<-FALSE)
	ifelse(!is.null(x$call$nsim)&&(x$call$nsim>0),ci<-TRUE,ci<-FALSE)
	if(is.null(x$call[["H0"]])||(x$call[["H0"]]=="pitor")) h0<-"PI-tor"
	else if(x$call[["H0"]]=="pimim") h0<-"PI-mim"
	else h0<-"RL"
	def.par <- par(no.readonly = TRUE)
	on.exit(par(def.par))
	#if(options()$device=="windows")
	#	csize<-0.75*csize
	opt<-opt[1]
	if(opt=="all")
		mylayout<-layout(matrix(c(1,1,1,1,2,2,3,3,2,2,3,3,4,4,5,5,4,4,5,5),ncol=4,byrow=TRUE))
	else if(opt%in%c("L","K","n","g"))
		mylayout<-layout(matrix(c(1,1,1,1,rep(2,16)),ncol=4,byrow=TRUE))
	else
		stopifnot(opt%in%c("all","L","K","n","g"))
	if(missing(cols))
		cols=c(1,2,3)
	else if(length(cols)!=3)
		cols=c(cols,cols,cols)
	if(missing(lty))
			lty=c(1,3,2)
	else if(length(lty)!=3)
		lty=c(lty,lty,lty)
	if(missing(main))
		main<-deparse(x$call,width.cutoff=100)		
	if(missing(sub))
		sub<-c("pair density function","second-order neighbour density function","intertype function","modified intertype function : sqrt[K12(r)/pi]-r")
	if(ci) {
		alpha<-x$call[["alpha"]]
		p<-ifelse(!is.null(alpha),signif(100*(1-alpha),digits=6),99)
		par(mar=c(0.1,0.1,0.1,0.1),cex=csize)
		plot(x$r,x$g12$obs/2,type="n",axes=FALSE,xlab="",ylab="")
		if(legend)
			legend("center",c("obs",paste("theo (",h0,")",sep=""),paste(p,"% CI of",h0)),cex=1.5,lty=lty[1:3],bty="n",horiz=TRUE,title=main,col=cols[1:3],...)
		else
			legend("center","",cex=1.5,bty="n",horiz=TRUE,title=main,...)
		par(mar=c(5,5,0.1,2),cex=ifelse(opt%in%c("all"),0.75*csize,csize))
		if(opt%in%c("all","g")) { # g12-function
			lim<-range(x$g12[,1:4])
			plot(x$r,x$g12$obs,ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[1]),type="n",xlab="distance (r)",ylab="g12(r)",cex.lab=1.25)
			lines(x$r,x$g12$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$g12$theo,lty=lty[2],col=cols[2],...)
			lines(x$r,x$g12$sup,lty=lty[3],col=cols[3],...)
			lines(x$r,x$g12$inf,lty=lty[3],col=cols[3],...)
		}
		if(opt%in%c("all","n")) { # n12-function
			lim<-range(x$n12[,1:4])
			plot(x$r,x$n12$obs,ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[2]),type="n",xlab="distance (r)",ylab="n12(r)",cex.lab=1.25)
			lines(x$r,x$n12$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$n12$theo,lty=lty[2],col=cols[2],...)
			lines(x$r,x$n12$sup,lty=lty[3],col=cols[3],...)
			lines(x$r,x$n12$inf,lty=lty[3],col=cols[3],...)
		}
		if(opt%in%c("all","K")) { # K-function
			plot(x$r,x$k12$obs,ylim=range(x$k12[,1:4]),main=paste("\n\n",sub[3]),type="n",xlab="distance (r)",ylab="K12(r)",cex.lab=1.25)
			lines(x$r,x$k12$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$k12$theo,lty=lty[2],col=cols[2],...)
			lines(x$r,x$k12$sup,lty=lty[3],col=cols[3],...)
			lines(x$r,x$k12$inf,lty=lty[3],col=cols[3],...)
		}
		if(opt%in%c("all","L")) { # L-function
			lim<-range(x$l12[,1:4])
			plot(x$r,x$l12$obs,ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[4]),type="n",xlab="distance (r)",ylab="L12(r)",cex.lab=1.25)
			lines(x$r,x$l12$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$l12$theo,lty=lty[2],col=cols[2],...)
			lines(x$r,x$l12$sup,lty=lty[3],col=cols[3],...)
			lines(x$r,x$l12$inf,lty=lty[3],col=cols[3],...)
		}
	}
	else {
		par(mar=c(0.1,0.1,0.1,0.1),cex=csize)
		plot(x$r,x$g12$obs/2,type="n",axes=FALSE,xlab="",ylab="")
		if(legend)
			legend("center",c("obs",paste("theo (",h0,")",sep="")),cex=1.5,lty=lty[1:2],bty="n",horiz=TRUE,title=main,col=cols[1:2],...)
		else
			legend("center","",cex=1.5,bty="n",horiz=TRUE,title=main,...)
		par(mar=c(5,5,0.1,2),cex=ifelse(opt%in%c("all"),0.75*csize,csize))
		if(opt%in%c("all","g")) { # g-function
			lim<-range(x$g12)
			plot(x$r,x$g12$obs,ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[1]),type="n",xlab="distance step (r)",ylab="g12(r)",cex.lab=1.25,...)
			lines(x$r,x$g12$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$g12$theo,lty=lty[2],col=cols[2],...)
		}
		if(opt%in%c("all","n")) { # n-function
			lim<-range(x$n12)
			plot(x$r,x$n12$obs,ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[2]),type="n",xlab="distance step (r)",ylab="n12(r)",cex.lab=1.25,...)
			lines(x$r,x$n12$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$n12$theo,lty=lty[2],col=cols[2],...)
		}
		if(opt%in%c("all","K")) { # k-function
			plot(x$r,x$k12$obs,ylim=range(x$k),main=paste("\n\n",sub[3]),type="n",xlab="distance step (r)",ylab="K12(r)",cex.lab=1.25,...)
			lines(x$r,x$k12$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$k12$theo,lty=lty[2],col=cols[2],...)
		}
		if(opt%in%c("all","L")) { # L-function
			lim<-range(x$l12)
			plot(x$r,x$l12$obs,ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[4]),type="n",xlab="distance step (r)",ylab="L12(r)",cex.lab=1.25,...)
			lines(x$r,x$l12$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$l12$theo,lty=lty[2],col=cols[2],...)
		}
	}	
}

plot.fads.kpqfun<-function (x,opt=c("L","K","n","g"),cols,lty,main,sub,legend=TRUE,csize=1,...) {
	na<-length(x$labpq)
	nf<-ceiling(sqrt(na))
	def.par <- par(no.readonly = TRUE)
	on.exit(par(def.par))
	#if(options()$device=="windows")
	#	csize<-0.75*csize
	mylayout<-layout(matrix(c(rep(1,nf),seq(2,((nf*nf)+1),1)),(nf+1),nf,byrow=TRUE))
	opt<-opt[1]
	if(opt=="g") {
		val<-x$gpq
		theo<-matrix(rep(rep(1,na),each=length(x$r)),ncol=na)
		ylab=paste("gpq(r)",sep="")
	}
	if(opt=="n") {
		val<-x$npq
		theo<-matrix(rep(rep(x$intensity,nf),each=length(x$r)),ncol=na)
		ylab=paste("npq(r)",sep="")
	}
	if(opt=="K") {
		val<-x$kpq
		theo<-matrix(rep(pi*x$r^2,na),ncol=na)
		ylab=paste("Kpq(r)",sep="")
	}
	if(opt=="L") {
		val<-x$lpq
		theo<-matrix(rep(rep(0,na),each=length(x$r)),ncol=na)
		ylab=paste("Lpq(r)",sep="")
	}
	if(missing(cols))
		cols=c(1,2)
	else if(length(cols)!=2)
		cols=c(cols,cols)
	if(missing(lty))
			lty=c(1,3)
	else if(length(lty)!=2)
		lty=c(lty,lty)
	if(missing(main))
		main<-deparse(x$call,width.cutoff=100)
	if(missing(sub))
		sub<-x$labpq
	lim<-range(val)
	par(mar=c(0.1,0.1,0.1,0.1),cex=csize)
	plot.default(val[,1],val[,2]/2,type="n",axes=FALSE,xlab="",ylab="")
	if(legend)
		legend("center",c("obs","theo (CSR/PI)"),cex=1.5,lty=lty[1:2],bty="n",horiz=TRUE,title=main,col=cols[1:2],...)
	else
			legend("center","",cex=1.5,bty="n",horiz=TRUE,title=main,...)
	par(mar=c(5,5,0.1,2),cex=0.66*csize)
	for(i in 1:na) {
		plot(x$r,val[,i],ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[i]),type="n",xlab="distance (r)",ylab=ylab,cex.lab=1.25,...)
		lines(x$r,val[,i],lty=lty[1],col=cols[1],...)
		lines(x$r,theo[,i],lty=lty[2],col=cols[2],...)
	}	
}

plot.fads.kp.fun<-function (x,opt=c("L","K","n","g"),cols,lty,main,sub,legend=TRUE,csize=1,...) {
	na<-length(x$labp)
	nf<-ceiling(sqrt(na))
	def.par <- par(no.readonly = TRUE)
	on.exit(par(def.par))
	#if(options()$device=="windows")
	#	csize<-0.75*csize
	mylayout<-layout(matrix(c(rep(1,nf),seq(2,((nf*nf)+1),1)),(nf+1),nf,byrow=TRUE))
	opt<-opt[1]
	if(opt=="g") {
		val<-x$gp.
		theo<-matrix(rep(rep(1,na),each=length(x$r)),ncol=na)
		ylab=paste("gp.(r)",sep="")
	}
	if(opt=="n") {
		val<-x$np.
		intensity<-sum(x$intensity)-x$intensity
		theo<-matrix(rep(intensity,each=length(x$r)),ncol=na)
		ylab=paste("np.(r)",sep="")
	}
	if(opt=="K") {
		val<-x$kp.
		theo<-matrix(rep(pi*x$r^2,na),ncol=na)
		ylab=paste("Kp.(r)",sep="")
	}
	if(opt=="L") {
		val<-x$lp.
		theo<-matrix(rep(rep(0,na),each=length(x$r)),ncol=na)
		ylab=paste("Lp.(r)",sep="")
	}
	if(missing(cols))
		cols=c(1,2)
	else if(length(cols)!=2)
		cols=c(cols,cols)
	if(missing(lty))
			lty=c(1,3)
	else if(length(lty)!=2)
		lty=c(lty,lty)
	if(missing(main))
		main<-deparse(x$call,width.cutoff=100)
	if(missing(sub))
		sub<-paste(x$labp,"-all others",sep="")
	lim<-range(val)
	par(mar=c(0.1,0.1,0.1,0.1),cex=csize)
	plot.default(val[,1],val[,2]/2,type="n",axes=FALSE,xlab="",ylab="")
	if(legend)
		legend("center",c("obs","theo (PI)"),cex=1.5,lty=lty[1:2],bty="n",horiz=TRUE,title=main,col=cols[1:2],...)
	else
		legend("center","",cex=1.5,bty="n",horiz=TRUE,title=main,...)
	par(mar=c(5,5,0.1,2),cex=0.66*csize)
	for(i in 1:na) {
		plot(x$r,val[,i],ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[i]),type="n",xlab="distance (r)",ylab=ylab,cex.lab=1.25,...)
		lines(x$r,val[,i],lty=lty[1],col=cols[1],...)
		lines(x$r,theo[,i],lty=lty[2],col=cols[2],...)
	}	
}

plot.fads.kmfun<-function (x,opt=c("all","K","g"),cols,lty,main,sub,legend=TRUE,csize=1,...) {
	ifelse(!is.null(x$call$nsim)&&(x$call$nsim>0),ci<-TRUE,ci<-FALSE)
	def.par <- par(no.readonly = TRUE)
	on.exit(par(def.par))
	#if(options()$device=="windows")
	#	csize<-0.75*csize
	opt<-opt[1]
	if(opt=="all")
		mylayout<-layout(matrix(c(1,1,1,1,rep(2,8),rep(3,8)),ncol=4,byrow=TRUE))
	else if(opt%in%c("K","g"))
		mylayout<-layout(matrix(c(1,1,1,1,rep(2,16)),ncol=4,byrow=TRUE))
	else
		stopifnot(opt%in%c("all","K","g"))
	if(missing(cols))
		cols=c(1,2,3)
	else if(length(cols)!=3)
		cols=c(cols,cols,cols)
	if(missing(lty))
			lty=c(1,3,2)
	else if(length(lty)!=3)
		lty=c(lty,lty,lty)
	if(missing(main))
		main<-deparse(x$call,width.cutoff=100)
	if(missing(sub))
		sub<-c("pair correlation function","mark correlation function")
	if(ci) {
		alpha<-x$call[["alpha"]]
		p<-ifelse(!is.null(alpha),signif(100*(1-alpha),digits=6),99)
		par(mar=c(0.1,0.1,0.1,0.1),cex=csize)
		plot(x$r,x$gm$obs/2,type="n",axes=FALSE,xlab="",ylab="")
		if(legend)
			legend("center",c("obs","theo (IM)",paste(p,"% CI of IM")),cex=1.5,lty=lty[1:3],bty="n",horiz=TRUE,title=main,col=cols[1:3],...)
		else
			legend("center","",cex=1.5,bty="n",horiz=TRUE,title=main,...)
		par(mar=c(5,5,0.1,2),cex=ifelse(opt%in%c("all"),0.85*csize,csize))
		if(opt%in%c("all","g")) { # gm-function
			lim<-range(x$gm[,1:4])
			plot(x$r,x$gm$obs,ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[1]),type="n",xlab="distance (r)",ylab="gm(r)",cex.lab=1.25,...)
			lines(x$r,x$gm$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$gm$theo,lty=lty[2],col=cols[2],...)
			lines(x$r,x$gm$sup,lty=lty[3],col=cols[3],...)
			lines(x$r,x$gm$inf,lty=lty[3],col=cols[3],...)	
		}
		if(opt%in%c("all","K")) { # K-function
			plot(x$r,x$km$obs,ylim=range(x$km[,1:4]),main=paste("\n\n",sub[2]),type="n",xlab="distance (r)",ylab="Km(r)",cex.lab=1.25,...)
			lines(x$r,x$km$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$km$theo,lty=lty[2],col=cols[2],...)
			lines(x$r,x$km$sup,lty=lty[3],col=cols[3],...)
			lines(x$r,x$km$inf,lty=lty[3],col=cols[3],...)
		}
	}
	else {
		par(mar=c(0.1,0.1,0.1,0.1),cex=csize)
		plot(x$r,x$gm$obs/2,type="n",axes=FALSE,xlab="",ylab="")
		if(legend)
			legend("center",c("obs","theo (IM)"),cex=1.5,lty=lty[1:2],bty="n",horiz=TRUE,title=main,col=cols[1:2],...)
		else
			legend("center","",cex=1.5,bty="n",horiz=TRUE,title=main,...)
		par(mar=c(5,5,0.1,2),cex=ifelse(opt%in%c("all"),0.85*csize,csize))
		if(opt%in%c("all","g")) { # g-function
			lim<-range(x$gm)
			plot(x$r,x$gm$obs,ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[1]),type="n",xlab="distance (r)",ylab="gm(r)",cex.lab=1.25,...)
			lines(x$r,x$gm$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$gm$theo,lty=lty[2],col=cols[2],...)
		}
		if(opt%in%c("all","K")) { # k-function
			plot(x$r,x$km$obs,ylim=range(x$km),main=paste("\n\n",sub[2]),type="n",xlab="distance (r)",ylab="Km(r)",cex.lab=1.25,...)
			lines(x$r,x$km$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$km$theo,lty=lty[2],col=cols[2],...)
		}
	}	
}

plot.fads.ksfun<-function (x,opt=c("all","K","g"),cols,lty,main,sub,legend=TRUE,csize=1,...) {
	ifelse(!is.null(x$call$nsim)&&(x$call$nsim>0),ci<-TRUE,ci<-FALSE)
	def.par <- par(no.readonly = TRUE)
	on.exit(par(def.par))
#if(options()$device=="windows")
#	csize<-0.75*csize
	opt<-opt[1]
	if(opt=="all")
	mylayout<-layout(matrix(c(1,1,1,1,rep(2,8),rep(3,8)),ncol=4,byrow=TRUE))
	else if(opt%in%c("K","g"))
	mylayout<-layout(matrix(c(1,1,1,1,rep(2,16)),ncol=4,byrow=TRUE))
	else
	stopifnot(opt%in%c("all","K","g"))
	
#	opt<-opt[1]
#	if(opt=="all")
#		mylayout<-layout(matrix(c(1,1,1,1,2,2,3,3,2,2,3,3,4,4,4,4,4,4,4,4),ncol=4,byrow=TRUE))
#	else if(opt%in%c("K","P","g"))
#		mylayout<-layout(matrix(c(1,1,1,1,rep(2,16)),ncol=4,byrow=TRUE))
#	else
#		stopifnot(opt%in%c("all","K","P","g"))
	if(missing(cols))
		cols=c(1,2,3)
	else if(length(cols)!=3)
		cols=c(cols,cols,cols)
	if(missing(lty))
		lty=c(1,3,2)
	else if(length(lty)!=3)
		lty=c(lty,lty,lty)
	if(missing(main))
		main<-deparse(x$call,width.cutoff=100)
	if(missing(sub))
		sub<-c("Standardized Shimatani non-cumulative (beta) function","Standardized Shimatani cumulative (alpha) function")
	if(ci) {
		alpha<-x$call[["alpha"]]
		p<-ifelse(!is.null(alpha),signif(100*(1-alpha),digits=6),99)
		par(mar=c(0.1,0.1,0.1,0.1),cex=csize)
		plot(x$r,x$gs$obs/2,type="n",axes=FALSE,xlab="",ylab="")
		if(legend)
			legend("center",c("obs","theo (RL)",paste(p,"% CI of RL",sep="")),cex=1.3,lty=lty[1:3],bty="n",horiz=TRUE,title=main,col=cols[1:3],...)
		else
			legend("center","",cex=1.5,bty="n",horiz=TRUE,title=main,...)
		par(mar=c(5,5,0.1,2),cex=ifelse(opt%in%c("all"),0.85*csize,csize))
		if(opt%in%c("all","g")) { # gs-function
			lim<-range(x$gs[,1:4])
			plot(x$r,x$gs$obs,ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[1]),type="n",xlab="distance (r)",ylab="gs(r)",cex.lab=1.25,...)
			lines(x$r,x$gs$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$gs$theo,lty=lty[2],col=cols[2],...)
			lines(x$r,x$gs$sup,lty=lty[3],col=cols[3],...)
			lines(x$r,x$gs$inf,lty=lty[3],col=cols[3],...)	
		}
		if(opt%in%c("all","K")) { # Ks-function
			plot(x$r,x$ks$obs,ylim=range(x$ks[,1:4]),main=paste("\n\n",sub[2]),type="n",xlab="distance (r)",ylab="Ks(r)",cex.lab=1.25,...)
			lines(x$r,x$ks$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$ks$theo,lty=lty[2],col=cols[2],...)
			lines(x$r,x$ks$sup,lty=lty[3],col=cols[3],...)
			lines(x$r,x$ks$inf,lty=lty[3],col=cols[3],...)
		}
	}
	else {
		par(mar=c(0.1,0.1,0.1,0.1),cex=csize)
		plot(x$r,x$gs$obs/2,type="n",axes=FALSE,xlab="",ylab="")
		if(legend)
		legend("center",c("obs","theo (RL)"),cex=1.3,lty=lty[1:2],bty="n",horiz=TRUE,title=main,col=cols[1:2],...)
		else
		legend("center","",cex=1.5,bty="n",horiz=TRUE,title=main,...)
		par(mar=c(5,5,0.1,2),cex=ifelse(opt%in%c("all"),0.85*csize,csize))
		if(opt%in%c("all","g")) { # gs-function
			lim<-range(x$gs)
			plot(x$r,x$gs$obs,ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[1]),type="n",xlab="distance (r)",ylab="gs(r)",cex.lab=1.25,...)
			lines(x$r,x$gs$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$gs$theo,lty=lty[2],col=cols[2],...)
		}
		if(opt%in%c("all","K")) { # ks-function
			plot(x$r,x$ks$obs,ylim=range(x$ks),main=paste("\n\n",sub[2]),type="n",xlab="distance (r)",ylab="Ks(r)",cex.lab=1.25,...)
			lines(x$r,x$ks$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$ks$theo,lty=lty[2],col=cols[2],...)
		}
	}	
}

plot.fads.krfun<-function (x,opt=c("all","K","g"),cols,lty,main,sub,legend=TRUE,csize=1,...) {
	ifelse(!is.null(x$call$nsim)&&(x$call$nsim>0),ci<-TRUE,ci<-FALSE)
	ifelse((is.null(x$call[["H0"]])||(x$call[["H0"]]=="rl")),h0<-"RL",h0<-"SE")
	def.par <- par(no.readonly = TRUE)
	on.exit(par(def.par))
#if(options()$device=="windows")
#	csize<-0.75*csize
	opt<-opt[1]
	if(opt%in%c("all"))
		mylayout<-layout(matrix(c(1,1,1,1,rep(2,8),rep(3,8)),ncol=4,byrow=TRUE))
	else if(opt%in%c("K","g"))
		mylayout<-layout(matrix(c(1,1,1,1,rep(2,16)),ncol=4,byrow=TRUE))
	else
		stopifnot(opt%in%c("all","K","g"))
	if(missing(cols))
		cols=c(1,2,3)
	else if(length(cols)!=3)
		cols=c(cols,cols,cols)
	if(missing(lty))
		lty=c(1,3,2)
	else if(length(lty)!=3)
		lty=c(lty,lty,lty)
	if(missing(main))
		main<-deparse(x$call,width.cutoff=100)
	if(missing(sub))
		sub<-c("Standardized Rao non-cumulative function","Standardized Rao cumulative function")
	if(ci) {
		alpha<-x$call[["alpha"]]
		p<-ifelse(!is.null(alpha),signif(100*(1-alpha),digits=6),99)
		par(mar=c(0.1,0.1,0.1,0.1),cex=csize)
		plot(x$r,x$gr$obs/2,type="n",axes=FALSE,xlab="",ylab="")
		if(legend)
			legend("center",c("obs",paste("theo (",h0,")",sep=""),paste(p,"% CI of",h0)),cex=1.3,lty=lty[1:3],bty="n",horiz=TRUE,title=main,col=cols[1:3],...)
		else
			legend("center","",cex=1.5,bty="n",horiz=TRUE,title=main,...)
		par(mar=c(5,5,0.1,2),cex=ifelse(opt%in%c("all"),0.85*csize,csize))
		if(opt%in%c("all","g")) { # gr-function
			lim<-range(x$gr[,1:4])
			plot(x$r,x$gr$obs,ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[1]),type="n",xlab="distance (r)",ylab="gr(r)",cex.lab=1.25,...)
			lines(x$r,x$gr$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$gr$theo,lty=lty[2],col=cols[2],...)
			lines(x$r,x$gr$sup,lty=lty[3],col=cols[3],...)
			lines(x$r,x$gr$inf,lty=lty[3],col=cols[3],...)	
		}
		if(opt%in%c("all","K")) { # Kr-function
			plot(x$r,x$kr$obs,ylim=range(x$kr[,1:4]),main=paste("\n\n",sub[2]),type="n",xlab="distance (r)",ylab="Kr(r)",cex.lab=1.25,...)
			lines(x$r,x$kr$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$kr$theo,lty=lty[2],col=cols[2],...)
			lines(x$r,x$kr$sup,lty=lty[3],col=cols[3],...)
			lines(x$r,x$kr$inf,lty=lty[3],col=cols[3],...)
		}
	}
	else {
		par(mar=c(0.1,0.1,0.1,0.1),cex=csize)
		plot(x$r,x$gr$obs/2,type="n",axes=FALSE,xlab="",ylab="")
		if(legend)
			legend("center",c("obs",paste("theo (",h0,")",sep="")),cex=1.3,lty=lty[1:3],bty="n",horiz=TRUE,title=main,col=cols[1:3],...)
		else
			legend("center","",cex=1.5,bty="n",horiz=TRUE,title=main,...)
		par(mar=c(5,5,0.1,2),cex=ifelse(opt%in%c("all"),0.85*csize,csize))
		if(opt%in%c("all","g")) { # gr-function
			lim<-range(x$gr)
			plot(x$r,x$gr$obs,ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[1]),type="n",xlab="distance (r)",ylab="gr(r)",cex.lab=1.25,...)
			lines(x$r,x$gr$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$gr$theo,lty=lty[2],col=cols[2],...)
		}
		if(opt%in%c("all","K")) { # kr-function
			plot(x$r,x$kr$obs,ylim=range(x$kr),main=paste("\n\n",sub[2]),type="n",xlab="distance (r)",ylab="Kr(r)",cex.lab=1.25,...)
			lines(x$r,x$kr$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$kr$theo,lty=lty[2],col=cols[2],...)
		}
	}	
}

plot.fads.kdfun<-function (x,opt=c("all","K","g"),cols,lty,main,sub,legend=TRUE,csize=1,...) {
	ifelse(!is.null(x$call$nsim)&&(x$call$nsim>0),ci<-TRUE,ci<-FALSE)
	def.par <- par(no.readonly = TRUE)
	on.exit(par(def.par))
#if(options()$device=="windows")
#	csize<-0.75*csize
	opt<-opt[1]
	if(opt%in%c("all"))
		mylayout<-layout(matrix(c(1,1,1,1,rep(2,8),rep(3,8)),ncol=4,byrow=TRUE))
	else if(opt%in%c("K","g"))
		mylayout<-layout(matrix(c(1,1,1,1,rep(2,16)),ncol=4,byrow=TRUE))
	else
		stopifnot(opt%in%c("all","K","g"))
	if(missing(cols))
		cols=c(1,2,3)
	else if(length(cols)!=3)
		cols=c(cols,cols,cols)
	if(missing(lty))
		lty=c(1,3,2)
	else if(length(lty)!=3)
		lty=c(lty,lty,lty)
	if(missing(main))
		main<-deparse(x$call,width.cutoff=100)
	if(missing(sub))
		sub<-c("Standardized Shen non-cumulative function","Standardized Shen cumulative function")
	if(ci) {
		alpha<-x$call[["alpha"]]
		p<-ifelse(!is.null(alpha),signif(100*(1-alpha),digits=6),99)
		par(mar=c(0.1,0.1,0.1,0.1),cex=csize)
		plot(x$r,x$gd$obs/2,type="n",axes=FALSE,xlab="",ylab="")
		if(legend)
			legend("center",c("obs","theo (SE)",paste(p,"% CI of SE")),cex=1.3,lty=lty[1:3],bty="n",horiz=TRUE,title=main,col=cols[1:3],...)
		else
			legend("center","",cex=1.5,bty="n",horiz=TRUE,title=main,...)
		par(mar=c(5,5,0.1,2),cex=ifelse(opt%in%c("all"),0.85*csize,csize))
		if(opt%in%c("all","g")) { # gd-function
			lim<-range(x$gd[,1:4])
			plot(x$r,x$gd$obs,ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[1]),type="n",xlab="distance (r)",ylab="gd(r)",cex.lab=1.25,...)
			lines(x$r,x$gd$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$gd$theo,lty=lty[2],col=cols[2],...)
			lines(x$r,x$gd$sup,lty=lty[3],col=cols[3],...)
			lines(x$r,x$gd$inf,lty=lty[3],col=cols[3],...)	
		}
		if(opt%in%c("all","K")) { # Kd-function
			plot(x$r,x$kd$obs,ylim=range(x$kd[,1:4]),main=paste("\n\n",sub[2]),type="n",xlab="distance (r)",ylab="Kd(r)",cex.lab=1.25,...)
			lines(x$r,x$kd$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$kd$theo,lty=lty[2],col=cols[2],...)
			lines(x$r,x$kd$sup,lty=lty[3],col=cols[3],...)
			lines(x$r,x$kd$inf,lty=lty[3],col=cols[3],...)
		}
	}
	else {
		par(mar=c(0.1,0.1,0.1,0.1),cex=csize)
		plot(x$r,x$gd$obs/2,type="n",axes=FALSE,xlab="",ylab="")
		if(legend)
			legend("center",c("obs","theo (SE)"),cex=1.3,lty=lty[1:3],bty="n",horiz=TRUE,title=main,col=cols[1:3],...)
		else
			legend("center","",cex=1.5,bty="n",horiz=TRUE,title=main,...)
		par(mar=c(5,5,0.1,2),cex=ifelse(opt%in%c("all"),0.85*csize,csize))
		if(opt%in%c("all","g")) { # gd-function
			lim<-range(x$gd)
			plot(x$r,x$gd$obs,ylim=c(lim[1],lim[2]+0.1*diff(lim)),main=paste("\n\n",sub[1]),type="n",xlab="distance (r)",ylab="gd(r)",cex.lab=1.25,...)
			lines(x$r,x$gd$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$gd$theo,lty=lty[2],col=cols[2],...)
		}
		if(opt%in%c("all","K")) { # Kd-function
			plot(x$r,x$kd$obs,ylim=range(x$kd),main=paste("\n\n",sub[2]),type="n",xlab="distance (r)",ylab="Kd(r)",cex.lab=1.25,...)
			lines(x$r,x$kd$obs,lty=lty[1],col=cols[1],...)
			lines(x$r,x$kd$theo,lty=lty[2],col=cols[2],...)
		}
	}	
}

plot.fads.mimetic<-function (x,opt=NULL,cols,lty,main,sub,legend=TRUE,csize=1,cex.main=1.5,pos,...) {
  if(missing(cols))
		cols=c(1,2)
	else if(length(cols)!=2)
		cols=c(cols,cols)
	if(missing(lty))
		lty=c(1,1)
	else if(length(lty)!=2)
		lty=c(lty,lty)
	if(missing(main)) {
		call<-match.call()
		main<-deparse(eval(eval(expression(call))[[2]][[2]])$call,width.cutoff=100)
	}
	plot(x$r,x$l$obs,ylim=range(rbind(x$l$obs,x$l$sim)),main=main,type="n",xlab="distance (r)",ylab="L(r)",cex.lab=1.25,cex.main=cex.main,...)
	lines(x$r,x$l$obs,lty=lty[1],col=cols[1],...)
	lines(x$r,x$l$sim,lty=lty[2],col=cols[2],...)
	if(legend)
		legend(pos,c("obs","sim"),cex=1.5,lty=lty[1:2],bty="n",horiz=TRUE,col=cols[1:2],...)
}


