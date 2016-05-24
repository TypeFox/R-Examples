seerHome="~/data/SEER"
pops=DF=NULL # to get rid of NOTES on package check about "no visible binding" 
load(file.path(seerHome,"73/pops.RData")) # this loads in pops
pym=NULL;pyf=NULL
for (i in 0:18) 
{   pym[i+1]=with(pops,sum(population[(popsex==1)&(popage==i)&(popyear<2009)]))
    pyf[i+1]=with(pops,sum(population[(popsex==2)&(popage==i)&(popyear<2009)])) }

load(file.path(seerHome,"73/lymyleuk.RData")) # this loads in DF
head(DF,2)
levels(DF$cancer)
table(DF$cancer)["CML"]
table(DF$ICD9)["2051"] # 100 extra not accounted for
table(DF$histo3)["9863"]+table(DF$histo3)["9875"]+table(DF$histo3)["9876"]+table(DF$histo3)["9945"]
table(DF$histo3)["9863"]+table(DF$histo3)["9875"]+table(DF$histo3)["9876"]
# d=DF[(DF$histo2==9863)&(DF$numprims==1)&(DF$yrdx<2009),] # paper used SEER 1973-2008
d=DF[(DF$cancer=="CML")&(DF$seqnum<2)&(DF$yrdx<2009),] # paper used SEER 1973-2008
d$agerec=as.numeric(cut(d$agedx,breaks=c(0,1,seq(5,85,5),130)))

m=hist(d$agerec[d$sex==1],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts
f=hist(d$agerec[d$sex==2],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts
age=c(0.5,3,seq(7.5,87.5,5))
datam=data.frame(age,cases=m,py=pym,incid=m/pym)[6:19,] 
dataf=data.frame(age,cases=f,py=pyf,incid=f/pyf)[6:19,]

if(length(grep("linux",R.Version()$os))) windows <- function( ... ) X11( ... )
if(length(grep("darwin",R.Version()$os))) windows <- function( ... ) quartz( ... )
windows(height=6,width=6)
par(mfrow=c(1,1),mar=c(4.7,6,3.3,1),oma=c(0,0,0,0),lwd=3,cex.lab=1.8,cex.axis=1.7,cex.main=1.8)
with(datam,plot(age,incid,log="y",xlab="Age",col='blue',pch=1,
                ylab="",cex=2,yaxt="n",ylim=c(1e-6,1e-4),xlim=c(22,87)))
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
title("U.S. CML Incidence 1973-2008")
legend(20,1.3e-4,c(paste("Males      k =",format(ssolm$par["k"], digits=2)),
                   paste( "Females  k =",format(ssolf$par["k"],digits=2))),
       col=c("blue","red"),pch=1:2,text.col=c("blue","red"),bty="n",cex=1.7)



