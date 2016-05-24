#FishRes--------------------------------2011-04-06
# Two linked populations (reserved and fished).
# Functions that define and run the model,
# with parameters shared globally for time series.
#-------------------------------------------------
local(envir=.PBSmodEnv,expr={
locale = sys.frame(sys.nframe() - 1) # local environment

# Make an expression (until part of PBSmodelling)
makeExpression <- function(s) { # s is a single string
	eval( parse(text=paste("expression(",s,")",collapse="",sep="")) ); };

# Harvest rate function (sinusoid across n periods)
hFun <- function(t,n,fmin,fmax) {
	fmin + ( 1 + sin(2*pi*(t/n)) )* ((fmax-fmin)/2); };

# Example density-dependent functions g(x)
g1 <- function(x,r,M,gam) {
	if (r<=0) return(rep(0,length(x))); 
	if (abs(gam)<1e-4) out <- (M/r)^x else out <- (1+((r/M)^gam-1)*x)^(-1/gam); 
	if (gam<=1e-4) {z <- 1/(1-(r/M)^gam); out[x>=z] <- 0; };
	return(out); };

# See g(x)
seegx <- function(){
	check(); # apply the standard checks
	getWinVal(scope="L"); act <- getWinAct()[1]; unpackList(pbs,scope="L");
	resetGraph(); expandGraph(mfrow=c(2,1),oma=c(0,0,0,0),mar=c(4,4,1,0.5),mgp=c(2.5,0.5,0),las=1);
	xmax <- xscale;
	x <- seq(0,xmax,len=100); n <- length(x); x1 <- 1;
	r <- switch(mod,rCont,rDisc); M <- switch(mod,Mcont,Mdisc);

	for (i in 1:2) {
		y <- g1(x,r,M,gam); y1 <- g1(x1,r,M,gam);
		if (i==2) y <- x*y;
		xlim <- range(x,na.rm=TRUE); ylim <- range(y,na.rm=TRUE);
		plot(0,0,type="n",xlim=xlim,ylim=ylim,xlab="x",ylab=makeExpression(switch(i,"g(x)","x*g(x)")),cex.lab=large);
		switch (i,abline(h=0,v=0,col="grey"),abline(a=0,b=M/r,col="red",lwd=3));
		lines(x,y,col="dodgerblue",lwd=3)
		lines(x=c(par()$usr[1],x1,x1),y=c(y1,y1,par()$usr[3]),col="darkgoldenrod",lty=5,lwd=2);
		points(x1,y1,col="gold",pch=15,cex=medium); points(x1,y1,col="blue",pch=0,cex=medium);
		#addLabel(.5,.95,paste("Model",i),col="salmon",cex=large); 
		addLabel(.925,switch(i,.925,.80),adj=1,paste("\\*g =",gam),vfont=c("sans serif","bold"),cex=large); };
};

# Model R1 Continuous for k = 0 (use ODE)
modR1a <- function(t,y,parms=NULL) {
	M <- Mcont; r <- rCont; a <- aCont;
	y1 <- y[1]; y2 <- y[2];
	F2 <- hFun(t,Fcy,FminCont,FmaxCont); F1 <- ifelse(t<tres,F2,FresCont);
	dy1 <- -M*y1 + r*y1*g1(y1/K1,r,M,gam) + a*((y2/p2) - (y1/p1)) - F1*y1;
	dy2 <- -M*y2 + r*y2*g1(y2/K2,r,M,gam) - a*((y2/p2) - (y1/p1)) - F2*y2;
	dy <- c(dy1=dy1,dy2=dy2);  # derivatives required for lsoda
	z <- c(dy1,dy2,F1,F2);        # values saved by lsoda
	names(z) <- c("dy1","dy2","F1","F2");
	list(dy,z)
}
# Model R1 Continuous for k > 0 (use DDE)
modR1b <- function(t,y,parms=NULL) {
	M <- Mcont; r <- rCont; a <- aCont;
	F2 <- hFun(t,Fcy,FminCont,FmaxCont); F1 <- ifelse(t<tres,F2,FresCont);
	lag <- if (t<=k) yinit else pastvalue(t-k);
	dy1  <- -M*y[1] + r*lag[1]*g1(lag[1]/K1,r,M,gam) + a*((y[2]/p2)-(y[1]/p1)) - F1*y[1];
	dy2  <- -M*y[2] + r*lag[2]*g1(lag[2]/K2,r,M,gam) - a*((y[2]/p2)-(y[1]/p1)) - F2*y[2];
	return(list(c(dy1,dy2), c(ttt=t,yy1=y[1],yy2=y[2],dy1=dy1,dy2=dy2,F1=F1,F2=F2)))
	#return(list(c(dy1,dy2), c(dy1=dy1,dy2=dy2,F1=F1,F2=F2)))
}
# Model R2 Discrete
modR2 <- function(tt,yy,parms=NULL) {
	#stop("model 2 still needs modification")
	M <- Mdisc; r <- rDisc; a <- aDisc
	N <- matrix(NA,nrow=length(tt),ncol=5,dimnames=list(1:length(tt),c("tt","y1","y2","F1","F2")))
	F2 <- hFun(tt,Fcy,FminDisc,FmaxDisc); F1 <- ifelse(tt<tres,F2,FresDisc);
	N[,"tt"]  <- tt; N[,"F1"] <- F1; N[,"F2"] <- F2;
	N[1,"y1"] <- yy[1]; N[1,"y2"] <- yy[2];
	for (i in 2:(tmax+1)) {                             # t in (k+1):tmax
		Nt1 <- N[i-1,"y1"]; Nt2 <- N[i-1,"y2"];          # i=t+1 here, so i-1=t
		if (i<=(k+1)) { Nk1 <- yy[1]; Nk2 <- yy[2]; }
		else { Nk1 <- N[i-1-k,"y1"]; Nk2 <- N[i-1-k,"y2"]; }      # i=t+1 here, so i-1=t
		N[i,"y1"] <- max(0, Nt1 -M*Nt1 + r*Nk1*g1(Nk1/K1,r,M,gam) + a*((Nt2/p2) - (Nt1/p1)) - N[i-1,"F1"]*Nt1);
		N[i,"y2"] <- max(0, Nt2 -M*Nt2 + r*Nk2*g1(Nk2/K2,r,M,gam) - a*((Nt2/p2) - (Nt1/p1)) - N[i-1,"F2"]*Nt2);
	}
	return(N)
}
check <- function() { 
	getWinVal(scope="L"); p1p2 <- p1*(1-p1);
	#convert model 1 parms to equivalent model 2 parms
	if (mod==1) { FresDisc <- 1-exp(-FresCont); FminDisc <- 1-exp(-FminCont); FmaxDisc <- 1-exp(-FmaxCont);
		hbig <- 1-exp(-Fbig); Mdisc <- 1-exp(-Mcont); rDisc <- exp(rCont)-1; aDisc <- p1p2*(1-exp(-aCont/p1p2));};
	if (mod==2) { FresCont = -log(GT0(1-FresDisc)); FminCont = -log(GT0(1-FminDisc))
		FmaxCont = -log(GT0(1-FmaxDisc)); Fbig = -log(GT0(1-hbig)); Mcont = -log(GT0(1-Mdisc))
		rCont = log(GT0(1+rDisc)); aCont = -p1p2*log(GT0(1-aDisc/p1p2)) }
	Fcont <- c(FresCont,FminCont,FmaxCont,Fbig);
	Fdisc <- c(FresDisc,FminDisc,FmaxDisc,hbig);
	
	if (enableCheck) {
		if (mod==1) {
			if (rCont <= 0) {rCont <- 0.1; cat("\n***** r reset to 0.1; reminder r > 0");};
			if (Mcont<=0 | Mcont>=rCont) {Mcont <- rCont/2; cat("\n***** M reset to mid-point; reminder 0 < M < r");};
			if (any(Fcont<0)) {Fcont<-pmax(Fcont,0,na.rm=TRUE); FresCont<-Fcont[1]; FminCont<-Fcont[2]; FmaxCont<-Fcont[3]; Fbig<-Fcont[4];}; 
			aCont <- pmax(0,aCont);
			FresDisc <- 1-exp(-FresCont); FminDisc <- 1-exp(-FminCont); FmaxDisc <- 1-exp(-FmaxCont);
			hbig <- 1-exp(-Fbig); Mdisc <- 1-exp(-Mcont); rDisc <- exp(rCont)-1; aDisc <- p1p2*(1-exp(-aCont/p1p2)); };
		if (mod==2) {
			k <- round(abs(k)); 
			if (rDisc <= 0) {rDisc <- 0.1; cat("\n***** r reset to 0.1; reminder r > 0");};
			if (Mdisc<=0 | Mdisc>=min(1,rDisc)) {Mdisc <- min(1,rDisc)/2; cat("\n***** M reset to mid-point; reminder 0 < M < min(1,r)");};
			aDisc <- pmax(0,aDisc);
			if (aDisc>=p1*(1-p1)) { aDisc <- p1*(1-p1)-.01; cat("\n***** a reset maximum-0.01; remember a < p1*p2"); };
			if (any(Fdisc<0)) {Fdisc<-pmax(Fdisc,0,na.rm=TRUE); FresDisc<-Fdisc[1]; FminDisc<-Fdisc[2]; FmaxDisc<-Fdisc[3]; hbig<-Fdisc[4];}; 
			if (any(Fdisc>1)) {Fdisc<-pmin(Fdisc,1,na.rm=TRUE); FresDisc<-Fdisc[1]; FminDisc<-Fdisc[2]; FmaxDisc<-Fdisc[3]; hbig<-Fdisc[4];};
			FresCont = -log(GT0(1-FresDisc)); FminCont = -log(GT0(1-FminDisc)); FmaxCont = -log(GT0(1-FmaxDisc)); 
			Fbig = -log(GT0(1-hbig)); Mcont = -log(GT0(1-Mdisc)); rCont = log(GT0(1+rDisc)); aCont = -p1p2*log(GT0(1-aDisc/p1p2)) }
	}
	alist <- list(FresCont=FresCont,FminCont=FminCont,FmaxCont=FmaxCont,
	              FresDisc=FresDisc,FminDisc=FminDisc,FmaxDisc=FmaxDisc,
	              Fbig=Fbig,hbig=hbig,rDisc=rDisc,rCont=rCont,Mdisc=Mdisc,
	              Mcont=Mcont,k=k, aDisc=aDisc, aCont=aCont);
	alist <- sapply(alist,round,3); setWinVal(alist)
}
#-----------------------------------------------------------------------
# Run Models R using shared global parameters
runModel <- function() {
	check();
	clrs <- c("forestgreen","red","blue","dodgerblue3");
	remove(list=ls(locale)[is.element(ls(locale),c(names(getWinVal()),"p2","K1","K2"))],pos=locale) # remove global values from previous runs
	getWinVal(scope="P"); #getWinVal(scope="L")
	act <- getWinAct()[1]; unpackList(pbs,scope="L")

	if(ptype=="p"){
		if (!exists("yout",envir=.PBSmodEnv)){
			resetGraph(); addLabel(.5,.5,"Run Time Series First",col="red",cex=large); return();}
		else {
			resetGraph();  par(mgp=c(2,.75,0),las=1) #; tget(yout)
			pairs(yout,pch=20,cex=tiny,col=clrs[4],gap=0,cex.labels=huge, cex.axis=big); return(); } };

	p2 <- 1-p1; K1 <- p1*K;  K2 <- (1-p1)*K;
	yinit <- c(p1*x10*K,(1-p1)*x20*K); # see Theorem 1
	pnams <- c("yinit","K1","K2","p2");  # extra parameters used in models, make global to speed up solvers
	plist <- list(); for (i in pnams) plist[[i]] <- get(i);
	unpackList(plist,scope="P");

	if (mod==1 && k==0) {
		tt <- seq(0,tmax,tstp);
		yout <- lsoda(y=yinit,times=tt,func=modR1a,parms=plist,rtol=rtol,atol=atol);
		colnames(yout) <- c("tt","y1","y2","dy1","dy2","F1","F2"); }
	if (mod==1 && k>0) {
		tt <- seq(0,tmax,tstp);
		yout <- dde(y=yinit,times=tt,func=modR1b,parms=plist,tol=atol); #tput(yout)
		colnames(yout)[1:3] <- c("tt","y1","y2"); }; 
	if (mod==2) {
		tt <- seq(0,tmax,by=1);
		yout <- modR2(tt,yinit,plist); };

	yout <- as.data.frame(yout); assign("yout",yout,pos=.PBSmodEnv);
	ptype <- getWinVal("ptype");
	unpackList(yout,scope="L"); xlim <- range(tt);

	resetGraph(); 
	if(ptype=="t") {
		expandGraph(mfrow=c(switch(mod,4,3),1),mar=c(0,5,0,0.5), oma=c(4,0,0.5,0), 
			mgp = c(3.25,0.3,0),las=1,tck=.02,cex.axis=big,cex.lab=huge);
		yfac <-.08;
		Nlim <- c(0,max(y1,y2,y1+y2)); Nlim <- Nlim + yfac*c(-1,1)*diff(Nlim);
		Flim <- range(F1,F2); Flim <- Flim + yfac*c(-1,1)*diff(Flim);
		C1 <- y1*F1; C2 <- y2*F2;
		C12 <- C1 + C2; Clim <- range(C1,C2,C12); Clim <- Clim + yfac*c(-1,1)*diff(Clim);

		plot(0,0,xlim=xlim,ylim=Nlim, xlab="",ylab="N",type="n");
		abline(h=K,col="grey");
		lines(tt,y1+y2,col=clrs[3],lwd=2);
		lines(tt,y1,col=clrs[1],lwd=2); lines(tt,y2,col=clrs[2],lwd=2);
		legend(x="topright",lwd=3,col=clrs[1:3],legend=c("Reserve","Fishery","Total"),cex=big,horiz=TRUE,bty="n")

		if (mod==1) {
			dNlim <- range(dy1,dy2,dy1+dy2); dNlim <- dNlim + yfac*c(-1,1)*diff(dNlim);
			plot(0,0,xlim=xlim,ylim=dNlim,xlab="",ylab="dN/dt",type="n");
			lines(tt,dy1+dy2,col=clrs[3],lwd=2);
			lines(tt,dy1,col=clrs[1],lwd=2); lines(tt,dy2,col=clrs[2],lwd=2); };

		plot(0,0,xlim=xlim,ylim=Flim,type="n",xlab="",ylab="F");
		lines(tt,F2,col=clrs[2],lwd=2); lines(tt,F1,col=clrs[1],lwd=2); 

		plot(0,0,xlim=xlim,ylim=Clim, xlab="",ylab="C",type="n");
		title("Cactch C",line=0.5,col.main="#400080");
		lines(tt,C12,col=clrs[3],lwd=2); lines(tt,C2,col=clrs[2],lwd=2); lines(tt,C1,col=clrs[1],lwd=2); 
		mtext("time",outer=TRUE,side=1,line=2.5,adj=.55,cex=big); }; 
	invisible(yout); };

#-------------------------------------------------------------
yield <- function (act=NULL) { # Equilibrium Yield Equations (Table 4)
	check();
	remove(list=ls(locale)[is.element(ls(locale),c("Yout","Zout"))], pos=locale)
	getWinVal(scope="L")
	if (is.null(act))  act <- getWinAct()[1]
	unpackList(pbs,scope="L")
	x1 <- seq(xvec[1],xvec[2],xvec[3]); nx <- length(x1);
	P1 <- seq(pvec[1],pvec[2],pvec[3]); np <- length(P1); 
	Fmax <- switch(mod,Fbig,hbig); hl <- switch(mod,c(FminCont,FmaxCont),c(FminDisc,FmaxDisc));
	a <- switch(mod, aCont, aDisc); 
	r <- switch(mod, rCont, rDisc);
	M <- switch(mod, Mcont, Mdisc);
	xlim<- c(0,1); ylim <- c(0,Fmax);
	off <- 0.02; xlim <- xlim + c(-1,1)*off*diff(xlim); ylim <- ylim + c(-1,1)*off*diff(ylim); 
	Yout <- array(NA,dim=c(nx,5,np),dimnames=list(1:nx,c("p1","F2","x1","x2","C2"),P1));
	clrs <- c("green3","orange","dodgerblue","forestgreen");
	
	for (p1 in P1) {
		pp <- as.character(p1);
		p2 <- 1-p1; K1 <- p1*K; K2 <- p2*K;
		Yout[,"p1",pp] <- rep(p1,nx); Yout[,"x1",pp] <- x1;
		#catch rate is the same for both models
		Q  <- (1 - (p1/a)*(r*g1(x1,r,M,gam)-M)); zQ <- Q<=0 | Q >1; Q[zQ] <- NA;
		x2 <- x1*Q; Yout[,"x2",pp] <- x2;
		F2 <- ((p1*x1)/(p2*x2)) * (r*g1(x1,r,M,gam)-M) + r*g1(x2,r,M,gam) - M; Yout[,"F2",pp] <- F2;
		C2 <- K2*F2*x2;  Yout[,"C2",pp] <- C2; #}
	};
	assign("Yout",Yout,pos=.PBSmodEnv); Zout <- NULL
	for (i in 1:np) {
		if (i==1) Zout <- Yout[,,i] else Zout <- rbind(Zout,Yout[,,i]); };
	Zout <- as.data.frame(Zout,row.names=1:nrow(Zout));
	bad <- apply(Zout,1,function(x){
		any(is.na(x)) | any(is.nan(x)) | any(is.infinite(x)); } );
	Zout <- Zout[!bad,];
	Zout <- Zout[Zout[,"F2"]>0 & Zout[,"F2"]<Fmax,];
	Zout <- Zout[rev(order(Zout$C2)),];
	setWinVal(list(eqOut=signif(Zout[1,c(5,1:4)],3)));
	assign("Zout",Zout,pos=.PBSmodEnv);
	CC <- Zout$C2; Cmax <- CC[1]; #CI <<- (1:length(CC))[rev(order(CC))];  ImaxC <<- CI[1];
	Z95 <- Zout[Zout$C2 >= (qlev * Cmax) & !is.na(Zout$C2),];

	if (any(act==c("contour","image"))) {
		zval <- c("C2","x1","x2");
		if (any(act==c("contour","image")))
			expandGraph(mfrow=c(3,1),mar=c(0,4.5,0,0.5), oma=c(3.5,0,0.5,0), mgp = c(2.75,0.3,0),
				cex.axis=big,cex.lab=huge,las=1,tck=.02);
		for (z in zval) {
			xx <- Zout$p1; yy <- Zout[,"F2"]; zz <- Zout[,z];
			x95 <- Z95$p1; y95 <- Z95[,"F2"];
			Zint <- interp(xx,yy,zz,duplicate="mean",linear=TRUE,
				xo=seq(min(xx),max(xx),length=ncell),yo=seq(min(yy),max(yy),length=ncell));
			#assign(paste("Zint",z,sep="."),Zint, pos=.PBSmodEnv);
			if (act=="contour")
				contour(Zint,nlevels=nlev,xlim=xlim,ylim=ylim,col=c(clrs[match(z,zval)],"grey50"),
					lwd=1,labcex=small,xlab="",ylab="F2");
			if (act=="image"){
				iclr <- rev(rainbow(40,start=0.05,end=0.7)); assign("zclr",iclr,pos=.PBSmodEnv)
				levs <- seq(ylim[1],ylim[2],len=length(iclr));
				zlev <- seq(min(zz,na.rm=TRUE),max(zz,na.rm=TRUE),len=length(iclr));
				ztck <- approx(zlev,levs,pretty(zz)); 
				ytck <- ztck$y[!is.na(ztck$y)]; ylab <- ztck$x[!is.na(ztck$y)]; ntck <-length(ytck);
				plot(0,0,type="n",xlim=xlim,ylim=ylim,xlab="",ylab="F2");
				image(Zint,add=TRUE,col=iclr,xlab="",ylab="");
				lines(x=c(par()$usr[1],pvec[2],NA,par()$usr[1],pvec[2]),y=c(hl[1],hl[1],NA,hl[2],hl[2]),col="grey10");
				points(x95,y95,pch=1,cex=medium);
				rect(1.02*xlim[2],levs[-length(levs)],1.04*xlim[2],levs[-1],col=iclr,border=FALSE);
				segments(rep(1.02*xlim[2],ntck),ytck,rep(1.04*xlim[2],ntck),ytck);
				text(1.01*xlim[2],ytck,ylab,adj=c(1,0),cex=medium);
			};
			points(xx[1],yy[1],col=ifelse(z=="C2" && act=="image","pink","red"),pch=17,cex=huge);
			points(xx[1],yy[1],col="black",pch=2,cex=huge); 
			xden <- apply(Zint$z,1,function(x){length(x[!is.na(x)])});
			if (xden[1]<xden[length(xden)]) space<-1 else space <-2;
			addLabel(switch(space,.07,.92),switch(space,.85,.07),z,cex=huge,col=clrs[match(z,zval)],adj=c(1,0));
			mtext("p1",outer=TRUE,side=1,line=2,cex=big,adj=.55); }; };
	
	if (act=="pairs") {
		resetGraph(); par(mgp=c(2,.75,0),las=1);
		clrs <- c("purple4","blue","green","yellow2","darkorange");
		CC <- Zout$C2;
		ccut <- as.numeric(cut(CC,quantile(CC,seq(0,1,.2))));  z <- clrs[ccut];
		pairs(Zout, gap=0, cex.labels=huge, cex.axis=big,
		panel=function(x,y,clr=z,...) {
			points(x,y,pch=16,cex=small,col=clr);
			points(x[1],y[1],col="red",pch=17,cex=large);
			points(x[1],y[1],col="black",pch=2,cex=large); } ); };
}; 

fig4 <- function(wmf=TRUE){
	check(); # apply the standard checks
	getWinVal(scope="L"); act <- getWinAct()[1]; unpackList(pbs,scope="L");
	resetGraph(); 
	if (wmf) win.metafile(filename="Fig4.wmf",width=8.5, height=10);
	expandGraph(mfrow=c(4,2),mar=c(2,5,0,0.5),oma=c(2,0,.5,0),mgp=c(3,0.6,0),las=1);
	xmax <- xscale;
	x <- seq(0,xmax,len=100); n <- length(x); x1 <- 1;
	r <- switch(mod,rCont,rDisc); M <- switch(mod,Mcont,Mdisc);
	gam <- c(-2,-1,0,2);
	
	for (j in gam) {
		for (i in 1:2) {
			y <- g1(x,r,M,j); y1 <- g1(x1,r,M,j);
			if (i==2) y <- x*y;
			xlim <- range(x,na.rm=TRUE); ylim <- range(y,na.rm=TRUE); ylim<-c(0,switch(i,1,0.42));
			plot(0,0,type="n",xlim=xlim,ylim=ylim,xlab="",ylab=makeExpression(switch(i,"g(x)","x*g(x)")),
				xaxt="n",cex.axis=big,cex.lab=huge);
			axis(1,labels=ifelse(par()$mfg[1]==4,TRUE,TRUE),cex.axis=big);
			if(par()$mfg[1]==4) mtext("x",side=1,line=2.5,cex=large);
			switch (i,abline(h=0,v=0,col="grey"),abline(a=0,b=M/r,col="red",lwd=3));
			lines(x,y,col="dodgerblue",lwd=3)
			lines(x=c(par()$usr[1],x1,x1),y=c(y1,y1,par()$usr[3]),col="darkgoldenrod",lty=5,lwd=2);
			points(x1,y1,col="gold",pch=15,cex=large); points(x1,y1,col="blue",pch=0,cex=large);
			#addLabel(.5,.95,paste("Model",i),col="salmon",cex=large); 
			addLabel(.925,switch(i,.85,.60),adj=1,paste("\\*g =",j),cex=huge,vfont=c("sans serif","plain")); }; };
	if (wmf) dev.off();
};

#===================================================================================================
pbs <- list(tiny=0.5, small=0.8, medium=1, big=1.25, large=1.5, huge=2); # cex sizes
remove(list=ls(locale)[is.element(ls(locale),c("yout","Yout","Zout"))], pos=locale) # remove objects from previous session
ips = installed.packages(); ip=ips[,"Package"]; names(ip)=ips[,"Version"]
for (i in c("PBSmodelling","deSolve","akima")) {
	if (any(ip==i)) { 
		eval(parse(text=paste("ipload=require(",i,", quietly=FALSE)",sep="")))
		if (!ipload) { 
			ii=paste("The",i,"package failed to load")
			showAlert(ii,i,"error"); stop(ii,call.=FALSE) } } 
	else {
		ii=paste("The `",i,"` package is required for this example.\n",sep="")
		if (getYes(paste(ii,"Load package from CRAN?",sep=""),"Package Needed")) install.packages(i)
		else stop(paste("Package `",i,"` needed for this example",sep=""),call.=FALSE)
	} }
		#showAlert(ii,i,"error"); stop(ii,call.=FALSE) } }
if (any(ip=="PBSddesolve")) 
	require(PBSddesolve,quietly=TRUE) else
	if (any(ip=="ddesolve")) {
		require(ddesolve,quietly=TRUE) 
		showAlert("ddesolve is now called PBSddesolve on CRAN","Information only") } else {
		ii="The PBSddesolve package is required for this example"
		showAlert(ii,"PBSddesolve","error"); stop(ii,call.=FALSE) } 
createWin("FishResWin.txt");

}) # end local scope
