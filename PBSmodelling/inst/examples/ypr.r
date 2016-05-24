# Ricker 1975 - p. 254 - Figure 10.2 (North Sea haddock from B-H)
# M=0.2, aref=1, Winf=1209, K=0.2, t0=-1.066, b=3
# Ricker 1975 - p. 257 - eqn. 10.24 (North Sea cod from Halliday)
# M=0.2, aref=t0, Winf=11.41, K=0.14, t0=0.07, b=3

local(envir=.PBSmodEnv,expr={
locale = sys.frame(sys.nframe() - 1) # local environment

# Incomplete beta function
incbeta <- function(x,pp,qq) {beta(pp,qq) * pbeta(x,pp,qq)}

ypr <- function(FF,tr) {
	#getWinVal(scope="L");
	r <- tr - t0; Z <- FF+M;
	y0 <- ((FF * Winf * exp(FF*r))/K) * incbeta( exp(-K*r), Z/K, b+1 );
	y0 * exp( M*(aref-t0) ); };

topt <- function(FF) { # optimal age for given F
	tr2 <- seq(tlim[1],tlim[2],length=tlim[3]*5); y2 <- ypr(FF,tr2);
	ym2 <- max(y2); tm2 <- mean(tr2[y2==ym2]);
	c(tm2,ym2); };

teum <- function(Flim) { #Eumetric line
	Fvec <- seq(Flim[1]+0.005,Flim[2],length=Flim[3]);
	Fmat <- matrix(Fvec,ncol=1);
	emat <- t(apply(Fmat,1,"topt"));
	zout <- cbind(Fvec,emat);
	dimnames(zout)[[2]] <- c("F","tR","Y");
	zout; };

RickFig <- function() {
	getWinVal(scope="P"); Fopt <- as.vector(Fin);
	Fvec <- seq(Flim[1],Flim[2],length=Flim[3]);
	tvec <- seq(tlim[1],tlim[2],length=tlim[3]);
	ymat <- outer(Fvec,tvec,FUN="ypr"); ymax <- max(ymat);
	EumP <- sapply(Fopt,topt); tput(EumP)
	Topt <- EumP[1,]; EumL <- rbind(Fopt,EumP);
	ye   <- teum(Flim);
	# clevs<- c(seq(clim[1],clim[2],by=clim[3]),floor(max(ye[,3])))
	clevs<- seq(clim[1],clim[2],by=clim[3]);
	clrs <- c("grey55","blue","red","red"); lty <- c(1,1,2,2);

	resetGraph(); expandGraph(mar=c(3,3,1,1),mgp=c(1.75,.5,0),las=1);
	contour(Fvec,tvec,ymat,levels=clevs,labcex=1.3,
		xlab="Fishing Mortality",ylab="Recruitment Age",lwd=2,col=clrs[1]);
	lines(ye[,1:2],lty=lty[2],lwd=2,col=clrs[2]);
	for (i in 1:length(Fopt)) {
		ii <- (i%%2+1)%%2+3;
		lines(c(Fopt[i],Fopt[i]),c(0,Topt[i]),lty=lty[ii],lwd=2,col=clrs[ii]);
		points(Fopt[i],Topt[i],pch=16,cex=1.5,col=clrs[ii]); 
		#text(Fopt[i],Topt[i],paste("(",paste(signif(EumL[,i],3),collapse=", "),")"),cex=.8,adj=c(-.1,1.1))
		};
	setWinVal(list(Rout=t(signif(EumP,5)))); };

autoC <- function() {
	getWinVal(scope="P");
	Fvec <- seq(Flim[1],Flim[2],length=Flim[3]);
	tvec <- seq(tlim[1],tlim[2],length=tlim[3]);
	ymat <- outer(Fvec,tvec,FUN="ypr"); yrng <- range(ymat,na.rm=T);
	nd   <- ifelse(yrng[2]<=1,2,ifelse(yrng[2]>1 & yrng[2]<=10,1,0));
	clim <- c(round(yrng,nd),round(diff(yrng)/nlevel,nd));
	setWinVal(list(clim=clim)); };	
	
#~~ <°)))<< ~~~~~~~~~~~~
require(PBSmodelling);
createWin("yprWin.txt");
}) # end local scope
