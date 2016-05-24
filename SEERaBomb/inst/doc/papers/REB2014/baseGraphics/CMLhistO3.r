# The point of this script is to discourage ICD-O3 9863 use because it is CML NOS (not otherwise specified).
# And to discourage O3 9875 (bcr-abl +) use because those tested may be selectively younger.
# Note that 9876 probably equals CMML since its aging rate constants are high. 
# Note also that ICD-O2 9863 is the sum of ICD-O3 9863 + 9875 + 9876

rm(list=ls(all=TRUE))
fexp<-function(p,dat)	{	c0=p[1];k=p[2]; 
	with(dat,{mn=exp(c0+k*age)*py
				-sum(cases*log(mn) - mn)})	}

load("/data/SEER/00/pops.RData") # this loads in pops
pym=NULL;pyf=NULL
for (i in 0:18) 
{   pym[i+1]=with(pops,sum(population[(popsex==1)&(popage==i)]))
	pyf[i+1]=with(pops,sum(population[(popsex==2)&(popage==i)])) }

load("/data/SEER/00/lymyleuk.RData") # this loads in DF
if(length(grep("linux",R.Version()$os))) windows <- function( ... ) X11( ... )
graphics.off()
windows(width=13,height=7,xpos=-300)
par(mfrow=c(2,3),mar=c(4.7,6,3.3,1),oma=c(0,0,0,0),lwd=3,cex.lab=1.8,cex.axis=1.7,cex.main=1.8)
for (i in 1:6) {
	if (i==1) {Indx=(DF$histo2==9863);tit="CML 9863 O2 "} 
	if (i==2) {Indx=(DF$histo2==9863)&(DF$seqnum<2);tit="CML 9863 O2 First"} 
	if (i==3) {Indx=(DF$histo2==9863)&(DF$numprims==1);tit="CML 9863 O2 Only 1"} 
	if (i==4) {Indx=(DF$histo3==9863);tit="CML NOS 9863 O3"} 
	if (i==5) {Indx=(DF$histo3==9863)&(DF$seqnum<2);tit="CML NOS 9863 O3 First"} 
	if (i==6) {Indx=(DF$histo3==9863)&(DF$numprims==1);tit="CML NOS 9863 O3 Only 1 "} 
	d=DF[Indx,] # paper used SEER 1973-2008
	m=hist(d$agerec[d$sex==1],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts
	f=hist(d$agerec[d$sex==2],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts
	age=c(0.5,3,seq(7.5,87.5,5))
	datam=data.frame(age,cases=m,py=pym,incid=m/pym)[6:19,] 
	dataf=data.frame(age,cases=f,py=pyf,incid=f/pyf)[6:19,]
	plot(age,m/pym,log="y",xlab="Age",type='p',col='blue',pch=1,ylab="",cex=2,yaxt="n",ylim=c(1e-7,1e-4)  )
	lines(age,f/pyf,type='p',col='red',cex=2,pch=2)  
	ssolm=optim(c(c0=-10,k=.1),fexp,dat=datam)
	ssolf=optim(c(c0=-10,k=.1),fexp,dat=dataf)
	ym=exp(ssolm$par["c0"]+ssolm$par["k"]*age[6:19])
	yf=exp(ssolf$par["c0"]+ssolf$par["k"]*age[6:19])
	axis(side=2,las=1, at=c(1e-7,1e-6,1e-5,1e-4),labels=expression(.1,1,10,10^2))
	lines(age[6:19],ym,col="blue")
	lines(age[6:19],yf,col="red")
	mtext(expression(paste("Cases per ",10^6," Person-Years")),side=2,line=3.3,cex=1.5)
	title(paste(tit,"2000-2010"))
	legend(-3,1.7e-4,
      c(sprintf("%3s Males; M/F=%3.2f; k = %5.3f",sum(m),mean(datam$incid/dataf$incid),ssolm$par["k"]),
		sprintf("%3s Females; k = %5.3f",sum(f),ssolf$par["k"])),
		col=c("blue","red"),pch=1:2,text.col=c("blue","red"),bty="n",cex=1.7)
}

windows(width=13,height=7,xpos=-300, ypos=-50)
par(mfrow=c(2,3),mar=c(4.7,6,3.3,1),oma=c(0,0,0,0),lwd=3,cex.lab=1.8,cex.axis=1.7,cex.main=1.8)
for (i in 1:6) {
	if (i==1) {Indx=(DF$histo3==9875);tit="CML 9875 O3 (BA-)"} 
	if (i==2) {Indx=(DF$histo3==9875)&(DF$seqnum<2);tit="CML 9875 O3 First"} 
	if (i==3) {Indx=(DF$histo3==9875)&(DF$numprims==1);tit="CML 9875 O3 Only 1 "} 
	if (i==4) {Indx=(DF$histo3==9876);tit="CML 9876 O3 (BA-)"} 
	if (i==5) {Indx=(DF$histo3==9876)&(DF$seqnum<2);tit="CML 9876 O3 First"} 
	if (i==6) {Indx=(DF$histo3==9876)&(DF$numprims==1);tit="CML 9876 O3 Only 1 "} 
	d=DF[Indx,] # paper used SEER 1973-2008
	m=hist(d$agerec[d$sex==1],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts
	f=hist(d$agerec[d$sex==2],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts
	age=c(0.5,3,seq(7.5,87.5,5))
	datam=data.frame(age,cases=m,py=pym,incid=m/pym)[6:19,] 
	dataf=data.frame(age,cases=f,py=pyf,incid=f/pyf)[6:19,]
	plot(age,m/pym,log="y",xlab="Age",type='p',col='blue',pch=1,ylab="",cex=2,yaxt="n",ylim=c(1e-7,1e-4)  )
	lines(age,f/pyf,type='p',col='red',cex=2,pch=2)  
	ssolm=optim(c(c0=-10,k=.1),fexp,dat=datam)
	ssolf=optim(c(c0=-10,k=.1),fexp,dat=dataf)
	ym=exp(ssolm$par["c0"]+ssolm$par["k"]*age[6:19])
	yf=exp(ssolf$par["c0"]+ssolf$par["k"]*age[6:19])
	axis(side=2,las=1, at=c(1e-7,1e-6,1e-5,1e-4),labels=expression(.1,1,10,10^2))
	lines(age[6:19],ym,col="blue")
	lines(age[6:19],yf,col="red")
	mtext(expression(paste("Cases per ",10^6," Person-Years")),side=2,line=3.3,cex=1.5)
	title(paste(tit,"2000-2010"))
	legend(-3,1.7e-4,
			c(sprintf("%3s Males; M/F=%3.2f; k = %5.3f",sum(m),mean(datam$incid/dataf$incid,na.rm=T),ssolm$par["k"]),
					sprintf("%3s Females; k = %5.3f",sum(f),ssolf$par["k"])),
			col=c("blue","red"),pch=1:2,text.col=c("blue","red"),bty="n",cex=1.7)
}
