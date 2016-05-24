# Figure 4 of current MS. CML incidence by race for 2000-2010
require(bbmle) 
load("/data/SEER/00/pops.RData") # this loads in pops
pyf=pym=vector(3,mode="list"); 
for (i in 0:18) { for (r in 1:2) {
		pym[[r]][i+1]=with(pops,sum(population[(popsex==1)&(popage==i)&(poprace==r)]))
		pyf[[r]][i+1]=with(pops,sum(population[(popsex==2)&(popage==i)&(poprace==r)]))}
	pym[[3]][i+1]=with(pops,sum(population[(popsex==1)&(popage==i)&(poprace>2)]))
	pyf[[3]][i+1]=with(pops,sum(population[(popsex==2)&(popage==i)&(poprace>2)])) }
p0=c(c0=-10,k=.04)
fexp<-function(p,dat){	c0=p[1];k=p[2]; 
	with(dat,{mn=exp(c0+k*(age[6:19]-55))*py[6:19]
				-sum(cases[6:19]*log(mn) - mn)})  }
nLL<-function(c0,k,x) with(x,{mn=exp(c0+k*(age-55))*py
				-sum(stats::dpois(cases, mn, log=TRUE))})

if(length(grep("linux",R.Version()$os))) windows <- function( ... ) X11( ... )
graphics.off();windows(width=12,height=5,xpos=-150)
par(mfrow=c(1,3),mar=c(4.7,0,0,0),oma=c(0,7.1,4,0.2),lwd=3,cex.lab=2.8,cex.axis=2.5,cex.main=2.8)
load("/data/SEER/00/lymyleuk.RData") # this loads in DF
for (i in 1:3) {
	if (i==3) d=DF[(DF$histo2==9863)&(DF$numprims==1)&(DF$race>2)&(DF$race<98),] else 
		d=DF[(DF$histo2==9863)&(DF$numprims==1)&(DF$race==i),] 
	N=sum(table(d$race))
	m=hist(d$agerec[d$sex==1],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts
	f=hist(d$agerec[d$sex==2],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts
	age=c(0.5,3,seq(7.5,87.5,5))
	plot(age,m/pym[[i]],log="y",xlab="Age",type='p',col='blue',pch=1,ylab="",cex=2,yaxt="n",ylim=c(2e-7,1.5e-4)  )
	lines(age,f/pyf[[i]],type='p',col='red',cex=2,pch=2)  
	ssolm=optim(p0,fexp,method="L-BFGS-B",dat=dm<-data.frame(age,cases=m,py=pym[[i]]))
	ssolf=optim(p0,fexp,method="L-BFGS-B",dat=df<-data.frame(age,cases=f,py=pyf[[i]]))
	ym=exp(ssolm$par["c0"]+ssolm$par["k"]*(age[6:19]-55))
	yf=exp(ssolf$par["c0"]+ssolf$par["k"]*(age[6:19]-55))
	if (i==1) axis(side=2,las=1, at=c(1e-6,1e-5,1e-4),labels=expression(1,10,10^2))
	lines(age[6:19],ym,col="blue"); lines(age[6:19],yf,col="red")
	MF=mean(ym/yf);       Tau=log(MF)/((ssolm$par["k"]+ssolf$par["k"])/2)  
	if (i==1) mtext(expression(paste("Cases per ",10^6," Person-Years")),side=2,line=3.8,cex=2)
	if (i==1) mtext(paste(N,"whites\nM/F =",format(MF,dig=2),"\nT =",format(Tau,dig=2),"Yrs"),side=1,line=-5.5,cex=1.5,adj=0.95)
	if (i==2) mtext(paste(N,"blacks\nM/F =",format(MF,dig=2),"\nT =",format(Tau,dig=2),"Yrs"),side=1,line=-5.5,cex=1.5,adj=0.95)
	if (i==3) mtext(paste(N,"others\nM/F =",format(MF,dig=2),"\nT =",format(Tau,dig=2),"Yrs"),side=1,line=-5.5,cex=1.5,adj=0.95)
	legend(-5,2.4e-4,c(paste("Males      k =",format(ssolm$par["k"],dig=2)),paste("Females  k =",format(ssolf$par["k"],dig=2))),
			col=c("blue","red"),pch=1:2,text.col=c("blue","red"),bty="n",cex=2) 
	fit0m <- mle2(nLL,start=list(c0=-10,k=.04),data=list(x=dm[6:19,]))
	fit0f <- mle2(nLL,start=list(c0=-10,k=.04),data=list(x=df[6:19,]))
	(CIF=cbind(coef(fit0f),confint(fit0f)))
	(CIM=cbind(coef(fit0m),confint(fit0m)))
	legend(43,2.4e-4,c(sprintf("(%s, %s)",format(CIM[2,2], digits=2),format(CIM[2,3], digits=2)),
					sprintf("(%s, %s)",format(CIF[2,2], digits=2),format(CIF[2,3], digits=2))),
			text.col=c("blue","red"),bty="n",cex=2)
}
title("CML Incidence in SEER 2000-2010",cex=3,outer=T)

