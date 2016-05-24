rm(list=ls(all=TRUE))
FIGS2=FALSE  # use 9863
FIGS2=TRUE   # use 205.1
load("/data/SEER/73/pops.RData") # this loads in pops
pops1=subset(pops,popyear<1985)
pops2=subset(pops,(popyear>=1985)&(popyear<1997))
pops3=subset(pops,popyear>=1997)
getSexPY<-function(pops) {pym=NULL;pyf=NULL
	attach(pops)
  for (i in 0:18) 
  {   pym[i+1]=sum(population[(popsex==1)&(popage==i)])
    pyf[i+1]=sum(population[(popsex==2)&(popage==i)]) }
  detach(pops)
  list(m=pym,f=pyf)
}
(py1=getSexPY(pops1))
(py2=getSexPY(pops2))
(py3=getSexPY(pops3))

load("/data/SEER/73/lymyleuk.RData") # this loads in DF
if (FIGS2) ndat=DF[(DF$ICD9==2051)&(DF$numprims==1),] else ndat=DF[(DF$histo3==9863)&(DF$numprims==1),] 
ndat1=subset(ndat,yrdx<1985)
ndat2=subset(ndat,(yrdx>=1985)&(yrdx<1997))
ndat3=subset(ndat,yrdx>=1997)

age=c(0.5,3,seq(7.5,87.5,5))
m1=with(ndat1,hist(agerec[sex==1],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts)
f1=with(ndat1,hist(agerec[sex==2],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts)
m2=with(ndat2,hist(agerec[sex==1],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts)
f2=with(ndat2,hist(agerec[sex==2],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts)
m3=with(ndat3,hist(agerec[sex==1],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts)
f3=with(ndat3,hist(agerec[sex==2],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts)

fexp<-function(p,dat){	c0=p[1];k=p[2]; 
	with(dat,{mn=exp(c0+k*age[6:19])*py[6:19]
				-sum(cases[6:19]*log(mn) - mn)})  }
p0=c(c0=-13,k=.04)  
ssolm=optim(p0,fexp,method="L-BFGS-B",control=list(maxit=400),dat=list(agex=age,cases=m1,py=py1$m))
ssolf=optim(p0,fexp,method="L-BFGS-B",control=list(maxit=400),dat=list(agex=age,cases=f1,py=py1$f))
print(km1<-ssolm$par["k"]); print(kf1<-ssolf$par["k"])  
ssolm=optim(p0,fexp,method="L-BFGS-B",control=list(maxit=400),dat=list(agex=age,cases=m2,py=py2$m))
ssolf=optim(p0,fexp,method="L-BFGS-B",control=list(maxit=400),dat=list(agex=age,cases=f2,py=py2$f))
print(km2<-ssolm$par["k"]); print(kf2<-ssolf$par["k"])  
ssolm=optim(p0,fexp,method="L-BFGS-B",control=list(maxit=400),dat=list(agex=age,cases=m3,py=py3$m))
ssolf=optim(p0,fexp,method="L-BFGS-B",control=list(maxit=400),dat=list(agex=age,cases=f3,py=py3$f))
print(km3<-ssolm$par["k"]); print(kf3<-ssolf$par["k"])  

graphics.off()
windows(height=5,width=12,xpos=-150)
par(mfrow=c(1,2),mar=c(4.7,6.2,0.8,1),lwd=3,cex.lab=2,cex.axis=2,cex.main=1.3)
plot(age,mincid1<-m1/py1$m,log="y",type='b',xlab="Age",pch=1,col='blue',ylab="",yaxt="n",ylim=c(0.9e-7,1.5e-4),cex=1.3)
lines(age,mincid2<-m2/py2$m,type='b',col='red',pch=2,cex=1.3)
lines(age,mincid3<-m3/py3$m,type='b',col='black',pch=3,cex=1.3)
axis(side=2,las=1, at=c(1e-7,1e-6,1e-5,1e-4),labels=expression(.1,1,10,10^2))
mtext(expression(paste("Cases per ",10^6," PY")),side=2,line=3.7,cex=2)
mtext("Males",side=3,line=-1.8,cex=2,adj=0.05)
legend(25,1.6e-6,c(paste("1973-1984   k =",format(km1,digits=2)), paste("1985-1996   k =",format(km2,digits=2)),
      paste("1997-2010   k =",format(km3,digits=2))), col=c("blue","red","black"),pch=1:3,bty="n",cex=1.5)

plot(age,fincid1<-f1/py1$f,log="y",type='b',xlab="Age",pch=1,cex=1.3,col='blue',ylab="",yaxt="n",ylim=c(.9e-7,1.5e-4)  )
lines(age,fincid2<-f2/py2$f,type='b',col='red',pch=2,cex=1.3)  
lines(age,fincid3<-f3/py3$f,type='b',col='black',pch=3,cex=1.3)  
mtext("Females",side=3,line=-1.8,cex=2,adj=0.05)
axis(side=2,las=1, at=c(1e-7,1e-6,1e-5,1e-4),labels=expression(.1,1,10,10^2))
mtext(expression(paste("Cases per ",10^6," PY")),side=2,line=3.5,cex=2)
legend(25,1.6e-6,c(paste("1973-1984   k =",format(kf1,digits=2)), paste("1985-1996   k =",format(kf2,digits=2)),
      paste("1997-2010   k =",format(kf3,digits=2))), col=c("blue","red","black"),pch=1:3,bty="n",cex=1.5)
