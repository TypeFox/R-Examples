#SGM------------------------------------2012-12-06
# Schnute Growth Model
#-------------------------------------------------
local(envir=.PBSmodEnv,expr={
locale = sys.frame(sys.nframe() - 1) # local environment

SGinit <- function(act=NULL) {
	SGlist <- readList("SGMdata.txt"); unpackList(SGlist,scope="P");
	if (is.null(act)) act <- getWinAct()[1]
	if (!is.null(act)) {
		pwr <- getWinVal()$pwr; 
		if (act=="reInit") {
			pwr <- 0; setWinVal(list(pwr=pwr)); };
		if (pwr>0) {
			mlen <- SGdata$len^pwr; SGdata$len <- mlen/max(mlen) * 100;
			mpars <- SGpars[3:4,1:3]^pwr; mpars<- round(mpars/max(mpars) * 100);
			mpars[1,1] <- max(mpars[1,1],1); SGpars[3:4,1:3] <- mpars
		}
	}
	tput(SGdata); tput(SGpars)
	setWinVal(list(parVec=SGpars)); };

SGfun <- function(P) {
	aModel  <- function(P,xobs,yobs,Pfix=NULL,is.pred=FALSE) {
		# Growth models - Schnute 1981
		a<-P[1]; b<-P[2]; y1<-P[3]; y2<-P[4];
		aa <- round(a,10); bb <- round(b,10); # For testing zero-values
		t1 <- Pfix[1]; t2 <- Pfix[2];
		#---Case 1---
		if (aa!=0 & bb!=0) {
			frac <- (1 - exp(-a*(xobs-t1))) / (1 - exp(-a*(t2-t1)));
			y0 <- y1^b + (y2^b - y1^b) * frac; y0 <- GT0(y0,eps=1e-8);
			y <- y0^(1/b); return(y);  }
		#---Case 2---
		if (aa!=0 & bb==0) {
			frac <- (1 - exp(-a*(xobs-t1))) / (1 - exp(-a*(t2-t1)));
			y <- y1 * exp(log(y2/y1) * frac); return(y);  }
		#---Case 3---
		if (aa==0 & bb!=0) {
			frac <- (xobs-t1) / (t2-t1);  y0 <- y1^b + (y2^b - y1^b) * frac;
			y0 <- GT0(y0,eps=1e-8); y <- y0^(1/b); return(y);  }
		#---Case 4---
		if (aa==0 & bb==0) {
			frac <- (xobs-t1) / (t2-t1);
			y <- y1 * exp(log(y2/y1) * frac); return(y);  }
	}
	tget(FP); unpackList(FP,scope="L");
	ypred <- aModel(P=P,xobs=xobs,yobs=yobs,Pfix=Pfix,is.pred=is.pred);
	if (is.pred) { FP$yobs <- ypred; tput(FP); return(ypred); }
	n <- length(yobs);  ssq <- sum( (yobs-ypred)^2 );
	return(n*log(ssq)); };

SGtest <- function(act=NULL) {
	getWinVal(scope="L"); tget(SGdata); tget(SGpars)
	FP <-list(Pfix=Pfix,is.pred=FALSE,xobs=SGdata$age,yobs=SGdata$len); tput(FP)
	Obag <- calcMin(pvec=parVec,func=SGfun,method=method,trace=trace,maxit=maxit,reltol=reltol,steptol=steptol,repN=repN); tput(Obag)
	tget(PBSmin); fmin <- PBSmin$fmin; np <- sum(parVec[,4]); ng <- nrow(SGdata);
	PBSmin$AICc <- 2*fmin + 2*np * (ng/(ng-np-1)); tput(PBSmin)
	P <- PBSmin$end; ftime <- PBSmin$time;
	Pcalc <- calcP(P,Pfix); tput(Pcalc) # t0, yinf, tstar, ystar, zstar (Eqns 24-28, Schnute 1981)
	Pcfig <- sapply(Pcalc,signif,5); Pfig <- signif(P,5);

	resetGraph(); expandGraph();
	clrs <- c("red","blue","green4","purple4","darkorange3","cornflowerblue","darkolivegreen");
	clr <- clrs[match(method,c("nlminb","nlm","Nelder-Mead","BFGS","CG","L-BFGS-B","SANN"))];
	xnew <- seq(SGdata$age[1],SGdata$age[ng],len=100);
	FP$is.pred <- TRUE;  FP$xobs <- xnew; tput(FP)
	ynew <- SGfun(P);
	plot(SGdata,las=1,cex.axis=1,mgp=c(2,0.5,0),xlab="Age",
		ylab=paste("Length",ifelse(pwr!=0 & pwr!=1,paste("^",round(pwr,3)),"")));
	axis(1,at=1:20,tck=-.01,labels=FALSE); axis(2,at=seq(0,100,5),tck=-.01,labels=FALSE);
	lines(xnew,ynew,col=clr,lwd=2);
	addLabel(.05,.95,paste("Method =",method),cex=1.2,adj=0,col=clr);
	addLabel(.05,.92,paste(paste(c("a","b","y1","y2"),Pfig,sep=" = "),collapse="\n"),adj=c(0,1),cex=0.8);
	addLabel(.22,.92,paste(paste(c("t0","yinf","t*","y*","z*"),Pcfig,sep=" = "),collapse="\n"),adj=c(0,1),cex=0.8);
	addLabel(.05,.80,paste("Timing =",paste(round(ftime[1],2),collapse=", "),"sec"),adj=0,cex=0.7,col="grey35");

	unpackList(Obag,scope="L"); unpackList(Pcalc,scope="L");
	Gbag <- list(Git=iters, Gev=evals, Gct=round(cpuTime,nd), Get=round(elapTime,nd), 
		Gf1=round(fminS,nd), Gf2=round(fminE,nd), Gaic=round(AIC,nd), Gaicc=round(PBSmin$AICc,nd),
		Gv1=round(Pend,nd), Gv2=sapply(Pcalc,round,nd), Gmess=message);
	setWinVal(Gbag); };

SGset <- function(){
	getWinVal(scope="L");
	pvec <- parVec; pvec[,1] <- Gv1;
	Gbag <- list(parVec=pvec,Git=0, Gev=0, Gct=0, Get=0, Gf1=0, Gf2=0, Gaic=0, Gaicc=0,
		Gv1=c(0,0,0,0), Gv2=c(0,0,0,0,0), 
		Gmess="------------\nparVec reset with last set of estimated parameters.\n------------");
	setWinVal(Gbag); };

calcP <- function(P,Pfix) { # t0, yinf, tstar, ystar, zstar (Eqns 24-28, Schnute 1981)
	a<-P[1]; b<-P[2]; y1<-P[3]; y2<-P[4]; t1 <- Pfix[1]; t2 <- Pfix[2];
	aa <- round(a,10); bb <- round(b,10);
	t0 <- NA; yinf <- NA; tstar <- NA; ystar <- NA;
	if (aa!=0 & bb!=0) {
		t0 <- t1 + t2 - (1/a)*log((exp(a*t2)*y2^b - exp(a*t1)*y1^b)/(y2^b - y1^b));
		yinf <- ((exp(a*t2)*y2^b - exp(a*t1)*y1^b)/(exp(a*t2)-exp(a*t1)))^(1/b);
		tstar <- t1 + t2 - (1/a)*log(b*(exp(a*t2)*y2^b - exp(a*t1)*y1^b)/(y2^b - y1^b));
		ystar <- ((1-b)*(exp(a*t2)*y2^b - exp(a*t1)*y1^b)/(exp(a*t2)-exp(a*t1)))^(1/b); };
	if (aa==0 & bb!=0) {
		t0 <- t1 + t2 - (t2*y2^b - t1*y1^b)/(y2^b - y1^b); };
	if (aa!=0 & bb==0) {
		yinf <- exp((exp(a*t2)*log(y2) - exp(a*t1)*log(y1))/(exp(a*t2)-exp(a*t1)));
		tstar <- t1 + t2 - (1/a)*log((exp(a*t2) - exp(a*t1))/log(y2/y1));
		ystar <- exp((exp(a*t2)*log(y2) - exp(a*t1)*log(y1))/(exp(a*t2)-exp(a*t1)) - 1); }
	zstar <- a / (1-b);
	return(list(t0=as.vector(t0),yinf=as.vector(yinf),tstar=as.vector(tstar),
		ystar=as.vector(ystar),zstar=as.vector(zstar))); };

show8 <- function(xy=1.5) {
	if (!exists("FP",where=.PBSmodEnv)) stop("Perform minimization first");
	tget(FP); unpackList(FP,scope="L"); P <- PBSmin$end;
	a<-P[1]; b<-P[2]; y1<-P[3]; y2<-P[4]; t1 <- Pfix[1]; t2 <- Pfix[2];
	bee  <- -a * (t2-t1) / log(y2/y1);
	xbee <- seq(-10,10,len=100); ybee <- bee*xbee;
	xone <- approx(x=ybee,y=xbee,xout=1)$y
	xlim <- c(-xy,xy); xlim[1] <- min(xlim[1],a,xone,na.rm=T); xlim[2] <- max(xlim[2],a,xone,na.rm=T);
	ylim <- c(-xy,xy); ylim[1] <- min(ylim[1],b,na.rm=T); ylim[2] <- max(ylim[2],b,na.rm=T);

	resetGraph(); expandGraph(mar=c(1,1,1,1));
	plot(0,0,type="n",xlab="",ylab="",xlim=xlim,ylim=ylim,axes=F);
	for (j in 1:2) { axis(j,pos=0,at=seq(-10,10,.5),las=1,cex.axis=.8,mgp=c(0,.6,0)); }
	abline(a=0,b=bee,lty=5,col="green4",lwd=2);
	segments(xone,1,xlim[2],1,lty=3,col="blue");
	points(a,b,pch=16,col="red",cex=1.5);
	text(a,b,paste("(",round(a,3),",",round(b,3),")",sep=""),cex=.9,adj=-.1); 
	addLabel(.025,.025,"See Fig. 1, p. 1133, Schnute (1981)\nCan. J. Fish. Aquat. Sci. 38:1128-1140",
		adj=c(0,0),cex=.7,col="grey30"); };

#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
if (!require(PBSmodelling, quietly=TRUE)) stop("The PBSmodelling package is required for this example")
createWin("SGMWin.txt")
remove(list=ls(locale)[is.element(ls(locale),c("FP","SGlist","SGdata","SGpars","PBSmin"))],pos=locale)
SGinit()
}) # end local scope
