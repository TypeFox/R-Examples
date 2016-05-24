###########################################
### Ecovirtual - Metapopulations Models ###
###########################################


### Propagulus Seed Rain 
metaPop <-function(cl,rw,fi,pc,pe, tmax)
{
	paisag=array(0,dim=c(rw,cl,tmax))
   nmanchas=cl*rw
	paisag[,,1]=matrix(sample(c(1,0),nmanchas,prob=c(fi,1-fi), replace=TRUE),rw,cl)
	resultado=numeric()
		for(tc in 2:tmax)
		{
	       paisag[,,tc][paisag[,,(tc-1)]==1]<-sample(c(0,1),sum(paisag[,,(tc-1)]), replace=TRUE, prob=c(pe,1-pe))
	       paisag[,,tc][paisag[,,(tc-1)]==0]<-sample(c(0,1),cl*rw-sum(paisag[,,(tc-1)]), replace=TRUE, prob=c(1-pc,pc))
	       resultado[tc-1]=sum(paisag[,,tc])/(cl*rw)	
	   }
          dev.new()
	animaMeta2(paisag)
	grFim(paisag)
	dev.new()
	F=pc/(pc+pe)
	plot(1:tmax,c(fi,resultado),type="l",xlab="Time",ylab="Proportion of occupation",
	ylim=c(0,1),main=paste("Propagulus rain","\n cols=",cl," rows=",rw," fi=",fi," pi=",pc," pe=",pe),font.lab=2,lwd=2)
	abline(h=F,col=2,lwd=2,lty=2)
	legend("topright", legend=("expected equilibrium"), lty=2, col="red", bty="n")
  invisible(paisag)
}

#metaPop(tmax=100,cl=20,rw=20,fi=0.2,pe=0.2,pc=0.5)


## Propagulus seed rain with Internal Colonization
metaCi <-function(cl,rw,fi,ci,pe, tmax)
{
paisag=array(0,dim=c(rw,cl,tmax))
nmanchas=cl*rw
paisag[,,1]=matrix(sample(c(rep(1,fi*nmanchas), rep(0,round((1-fi)*nmanchas)))), ncol=cl)
resultado=numeric()
	for(tc in 2:tmax)
	{
	pc=ci*sum(paisag[,,tc-1])/(cl*rw)
	paisag[,,tc][paisag[,,(tc-1)]==1]<-sample(c(0,1),sum(paisag[,,(tc-1)]), replace=TRUE,prob=c(pe,1-pe))
	paisag[,,tc][paisag[,,(tc-1)]==0]<-sample(c(0,1),cl*rw-sum(paisag[,,(tc-1)]), replace=TRUE,prob=c(1-pc,pc))
   resultado[tc-1]=sum(paisag[,,tc])/nmanchas
   }
dev.new()
animaMeta2(paisag)
grFim(paisag)
dev.new()
F=1-(pe/ci)
plot(1:tmax,c(fi,resultado),type="l",xlab="Time",ylab="Proportion of occupation",
ylim=c(0,1),main=paste("Propagulus Rain and Internal Colonization","\n cols=",cl," rows=",rw," fi=",fi," ci=",ci," pe=",pe),font.lab=2,lwd=2)
abline(h=F,col=2,lwd=2,lty=2)
legend("topright", legend=("expected equilibrium"), lty=2, col="red", bty="n")
invisible(paisag)
}

#metaCi(tmax=100,cl=10,rw=10,fi=.1,ci=1,pe=0.5)


## Propagulus Seed Rain with Rescue EFfect
metaEr <-function(cl,rw,fi,pc,ce, tmax)
{
nmanchas=cl*rw
paisag=array(0,dim=c(rw,cl,tmax))
paisag[,,1]=matrix(sample(c(1,0),nmanchas,prob=c(fi,1-fi), replace=TRUE),rw,cl)
resultado=numeric()
res=numeric()
	for(tc in 2:tmax)
	{
	pe=ce*(1-sum(paisag[,,tc-1])/nmanchas)
	paisag[,,tc][paisag[,,(tc-1)]==1]<-sample(c(0,1),sum(paisag[,,(tc-1)]), replace=TRUE, prob=c(pe,1-pe))
	paisag[,,tc][paisag[,,(tc-1)]==0]<-sample(c(0,1),cl*rw-sum(paisag[,,(tc-1)]), replace=TRUE, prob=c(1-pc,pc))
	resultado[tc-1]=sum(paisag[,,tc])/nmanchas
	res[tc-1]=pe
	}
dev.new()
animaMeta2(paisag)
grFim(paisag)
dev.new()
F=pc/ce
if(F>1){F=1}
pe.eq=ce-pc
if(pe.eq<0){pe.eq=0}
plot(1:tmax,c(fi,resultado),type="l",xlab="Time",ylab="Proportion/Probability", ylim=c(0,1),main=paste("Propagulus Rain and Rescue Effect","\n cols=",cl," rows=",rw," fi=",fi," pc=",pc," ce=",ce),font.lab=2,lwd=2) 
abline(h=F,col=2,lwd=2,lty=2) # equilibrio F
points(1:tmax,c(ce*(1-fi),res),type='l',lwd=2,col="blue") # pe observado
abline(h=pe.eq,col="green",lwd=2,lty=2) # pe equilibrio

ymin=min(resultado[(length(resultado)/2):(length(resultado))])
legend(x=length(resultado)/2,y=ymin, legend=c("Proportion of occupancy (P)", "Equilibrium P", "Extintion probability (pe)", "pe equilibrium"), lty=c(1,2,1,2), col=c("black","red","blue", "green"), bty="n")
invisible(paisag)
}

#metaEr(100,20,20,0.25,0.1,0.1)

## Propagulus Seed Rain with Internal Colonization and Rescue Effect
metaCiEr <-function(cl,rw,fi,ci,ce, tmax)
{
nmanchas=cl*rw
paisag=array(0,dim=c(rw,cl,tmax))
paisag[,,1]=sample(c(rep(0,round(nmanchas-fi*nmanchas)),rep(1,round(fi*nmanchas))))
resultado=numeric()
rese=numeric()
resi=numeric()
	for(tc in 2:tmax)
	{
	pe=ce*(1-(sum(paisag[,,tc-1])/nmanchas))
	pc=ci*sum(paisag[,,tc-1])/nmanchas
	paisag[,,tc][paisag[,,(tc-1)]==1]<-sample(c(0,1),sum(paisag[,,tc-1]),replace=TRUE,prob=c(pe,1-pe))
	paisag[,,tc][paisag[,,(tc-1)]==0]<-sample(c(0,1),nmanchas-sum(paisag[,,tc-1]),replace=TRUE,prob=c(1-pc,pc))
	resultado[tc-1]=sum(paisag[,,tc])/nmanchas
	rese[tc-1]=pe
	resi[tc-1]=pc
	}
dev.new()
animaMeta2(paisag)
grFim(paisag)
dev.new()
plot(1:tmax,c(fi,resultado),type="l",xlab="Time",ylab="Occupancy proportion", ylim=c(0,1),main=paste("Propagulus Rain and Internal colonization and Rescue Effect","\n cols=",cl," rows=",rw," fi=",fi," ci=",ci, "ce=",ce),font.lab=2,lwd=2)
abline(h=0,lty=2)
points(1:tmax,c(ce*(1-fi),rese),type='l',lwd=2,col=4,lty=3)
points(1:tmax,c(ci*fi,resi),type='l',lwd=2,col=6,lty=3)
legend("topright", legend=c("patchs occupancy", "colonization", "extinction"), lty=c(1,3,3), col=c(1,6,4), bty="n")
invisible(paisag)
}

#metaCiEr(100,10,10,0.5,0.5,0.5)
