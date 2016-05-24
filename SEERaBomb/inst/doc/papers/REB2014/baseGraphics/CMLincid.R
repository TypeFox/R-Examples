# Figure 1 of current MS. Differs from Blood Fig 5a in using 2000-2010, not 1973-2008
seerHome="/data/SEER"
load(file.path(seerHome,"00/pops.RData")) # this loads in pops
pym=NULL;pyf=NULL
for (i in 0:18) 
{ pym[i+1]=with(pops,sum(population[(popsex==1)&(popage==i)]))
	pyf[i+1]=with(pops,sum(population[(popsex==2)&(popage==i)])) }

load(file.path(seerHome,"00/lymyleuk.RData")) # this loads in DF
d=DF[(DF$histo2==9863)&(DF$numprims==1),] 
m=hist(d$agerec[d$sex==1],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts
f=hist(d$agerec[d$sex==2],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts
age=c(0.5,3,seq(7.5,87.5,5))
datam=data.frame(age,cases=m,py=pym,incid=m/pym)[6:19,] 
dataf=data.frame(age,cases=f,py=pyf,incid=f/pyf)[6:19,]

if(length(grep("linux",R.Version()$os))) windows <- function( ... ) X11( ... )
graphics.off();windows(width=6,height=6,xpos=-150)
par(mfrow=c(1,1),mar=c(4.7,6,3.3,1),oma=c(0,0,0,0),lwd=3,cex.lab=1.8,cex.axis=1.7,cex.main=1.8)
with(datam,plot(age,incid,log="y",xlab="Age",col='blue',pch=1,
				ylab="",cex=2,yaxt="n",ylim=c(2e-6,1e-4),xlim=c(22,87)))
with(dataf,points(age,incid,col='red',cex=2,pch=2))

fexp<-function(p,dat)	{	c0=p[1];k=p[2]; 
	with(dat,{mn=exp(c0+k*age)*py
				-sum(cases*log(mn) - mn)})	}
ssolm=optim(c(c0=-10,k=.04),fexp,method="L-BFGS-B",dat=datam)
ssolf=optim(c(c0=-10,k=.04),fexp,method="L-BFGS-B",dat=dataf)
ym=exp(ssolm$par["c0"]+ssolm$par["k"]*age[6:19])
yf=exp(ssolf$par["c0"]+ssolf$par["k"]*age[6:19])
axis(side=2,las=1, at=c(1e-6,1e-5,1e-4),labels=expression(1,10,10^2))
lines(age[6:19],ym,col="blue")
lines(age[6:19],yf,col="red")
mtext(expression(paste("Cases per ",10^6," Person-Years")),side=2,line=3.5,cex=2)
title("CML Incidence: SEER 2000-2010")
legend(19,1.3e-4,c(paste("Males      k =",format(ssolm$par["k"], digits=4)),
				paste( "Females  k =",format(ssolf$par["k"],digits=4))),
		col=c("blue","red"),pch=1:2,text.col=c("blue","red"),bty="n",cex=1.4)

# double check with mle2 and get k CI
require(bbmle) 
nLL<-function(c0,k,x) with(x,{mn=exp(c0+k*age)*py
				-sum(stats::dpois(cases, mn, log=TRUE))})
fit0m <- mle2(nLL,start=list(c0=-10,k=.04),data=list(x=datam))
fit0f <- mle2(nLL,start=list(c0=-10,k=.04),data=list(x=dataf))
(CIF=cbind(coef(fit0f),confint(fit0f)))
(CIM=cbind(coef(fit0m),confint(fit0m)))

legend(55,1.3e-4,c(sprintf("(%s, %s)",format(CIM[2,2], digits=3),format(CIM[2,3], digits=3)),
				sprintf("(%s, %s)",format(CIF[2,2], digits=3),format(CIF[2,3], digits=3))),
		col=c("blue","red"),text.col=c("blue","red"),bty="n",cex=1.4)

