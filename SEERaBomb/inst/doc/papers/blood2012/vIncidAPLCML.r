# This reproduces Fig S4 in  Radivoyevitch et al Blood 119, 4363-4371 (2012)
rm(list=ls(all=TRUE))
load("~/data/seer/73/pops.RData") # this loads in pops
pym=NULL;pyf=NULL
for (i in 0:18) 
{   pym[i+1]=with(pops,sum(population[(popsex==1)&(popage==i)&(popyear<2009)]))
	pyf[i+1]=with(pops,sum(population[(popsex==2)&(popage==i)&(popyear<2009)])) }

mylog<-function(x,fitEXP) if (fitEXP) x else log(x) 
fexp<-function(p,dat,fitEXP)	{	c0=p[1];k=p[2]; 
	with(dat,{mn=exp(c0+k*mylog(age,fitEXP))*py 
				-sum(cases*log(mn) - mn)})	}

load("~/data/SEER/73/lymyleuk.RData") # this loads in DF
# morphCodes=c(CML=9863,APL=9866)
if(length(grep("linux",R.Version()$os))) windows <- function( ... ) X11( ... )
if(length(grep("darwin",R.Version()$os))) windows <- function( ... ) quartz( ... )
graphics.off();
windows(width=8,height=8)
par(mfrow=c(2,2),mar=c(4.7,0,2.3,0),oma=c(4.4,6,0.4,0.4),lwd=3,cex.lab=1.8,cex.axis=1.7,cex.main=1.8)
for (fitEXP in c(TRUE,FALSE))
# 	for (i in names(morphCodes)) {
	for (i in c("CML","APL")) {
		d=DF[(DF$cancer==i)&(DF$seqnum<2)&(DF$yrdx<2009),]
		d$agerec=as.numeric(cut(d$agedx,breaks=c(0,1,seq(5,85,5),130)))
		m=hist(d$agerec[d$sex==1],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts
		f=hist(d$agerec[d$sex==2],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts
		age=c(0.5,3,seq(7.5,87.5,5))
		plot(age,m/pym,log="y",xlab="Age",col='blue',ylim=c(5e-7,1e-4),yaxt="n",ylab="",main=i)
		points(age,f/pyf,col='red'); 
		if (i=="CML") {
			mtext(expression(paste("Cases per ",10^6," Person-Years")),side=2,line=3.6,cex=1.5)
			axis(side=2,las=1,at=c(5e-7,1e-6,5e-6,1e-5,5e-5),labels=c(0.5,1,5,10,50))}
		datam=data.frame(age,cases=m,py=pym,incid=m/pym)[6:19,] 
		dataf=data.frame(age,cases=f,py=pyf,incid=f/pyf)[6:19,]
		p0=c(c0=-10,k=ifelse(fitEXP,.04,1))
# 		ssolm=optim(p0,fexp,method="L-BFGS-B",dat=datam,fitEXP=fitEXP)
# 		ssolf=optim(p0,fexp,method="L-BFGS-B",dat=dataf,fitEXP=fitEXP)
		ssolm=optim(p0,fexp,dat=datam,fitEXP=fitEXP)
		ssolf=optim(p0,fexp,dat=dataf,fitEXP=fitEXP)
		ym=exp(ssolm$par["c0"]+ssolm$par["k"]*mylog(age[6:19],fitEXP))
		yf=exp(ssolf$par["c0"]+ssolf$par["k"]*mylog(age[6:19],fitEXP))
		lines(age[6:19],ym,col="blue"); lines(age[6:19],yf,col="red")
		km<-format(ssolm$par["k"], digits=2);kf<-format(ssolf$par["k"], digits=2)
		mtext(ifelse(fitEXP,expression(y==Ae^{k%.%age}),expression(y==A%.%age^k)),side=3,line=-2.5,cex=1.3,adj=.05)
		mtext(bquote(k[m] == .(km)),side=3,line=-4,cex=1.3,adj=.05,col="blue")
		mtext(bquote(k[f] == .(kf)),side=3,line=-5.7,cex=1.3,adj=.05,col="red")
	} # the i loop
legend("bottomright",c("Males","Females"),col=c("blue","red"),cex=1.2,pch=1,bty="n")



