local(envir=.PBSmodEnv,expr={
locale = sys.frame(sys.nframe() - 1) # local environment

VBinit <- function() {
	VBlist <- readList("vonBdata.txt"); unpackList(VBlist,scope="P");
	setWinVal(list(parVec=VBpars)); };

VBfun <- function(P) {
	Linf <- P[1]; K <- P[2]; t0 <- P[3];
	obs  <- VBdata$len;
	pred <- Linf * (1 - exp(-K*(VBdata$age-t0)) );
	n    <- length(obs);
	ssq  <- sum( (obs-pred)^2 );
	return(n*log(ssq)); };

VBtest <- function() {
	getWinVal(scope="L");
	Obag <- calcMin(pvec=parVec,func=VBfun,method=method,trace=trace,maxit=maxit,reltol=reltol,steptol=steptol,repN=repN);
	tget(PBSmin); fmin <- PBSmin$fmin; np <- sum(parVec[,4]); ng <- nrow(VBdata);
	PBSmin$AICc <- 2*fmin + 2*np * (ng/(ng-np-1)); tput(PBSmin);
	P <- PBSmin$end; ftime <- PBSmin$time;

	resetGraph(); expandGraph();
	clrs <- c("red","blue","green4","purple4","darkorange3","cornflowerblue","darkolivegreen");
	clr <- clrs[match(method,c("nlminb","nlm","Nelder-Mead","BFGS","CG","L-BFGS-B","SANN"))];
	xnew <- seq(VBdata$age[1],VBdata$age[ng],len=100);
	ynew <- P[1] * (1 - exp(-P[2]*(xnew-P[3])) )
	plot(VBdata); lines(xnew,ynew,col=clr,lwd=2); 
	addLabel(.05,.95,paste("Method =",method),cex=1.2,adj=0,col=clr);
	addLabel(.05,.88,paste(paste(c("Linf","K","t0"),round(P,c(2,4,4)),sep=" = "),collapse="\n"),adj=0,cex=0.9);
	addLabel(.05,.80,paste("Timing =",paste(round(ftime[1],2),collapse=", "),"sec"),adj=0,cex=0.7,col="grey35");

	unpackList(Obag,scope="L");
	Gbag <- list(Git=iters, Gev=evals, Gct=round(cpuTime,nd), Get=round(elapTime,nd), 
		Gf1=round(fminS,nd), Gf2=round(fminE,nd), Gaic=round(AIC,nd), Gaicc=round(PBSmin$AICc,nd),
		Gv1=Pstart, Gv2=round(Pend,nd), Gmess=message);
	setWinVal(Gbag); };

VBset <- function(){
	getWinVal(scope="L");
	pvec <- parVec; pvec[,1] <- Gv2;
	Gbag <- list(parVec=pvec,Git=0, Gev=0, Gct=0, Get=0, Gf1=0, Gf2=0, Gaic=0, Gaicc=0,
		Gv1=c(0,0,0), Gv2=c(0,0,0), 
		Gmess="------------\nparVec reset with last set of estimated parameters.\n------------");
	setWinVal(Gbag); };

#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

require("PBSmodelling");
createWin("vonBWin.txt");
remove(list=ls(locale)[is.element(ls(locale),c("VBlist","VBdata","VBpars","PBSmin"))],envir=locale); VBinit()
}) # end local scope
