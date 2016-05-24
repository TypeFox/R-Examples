#Example where the juvenile stage duration is dependent upon temperature.

library(stagePop)

#solver.options=list(DDEsolver='deSolve',atol=1e-6,rtol=1e-6,hbsize=1e5)
solver.options=list(DDEsolver='PBS',tol=1e-8,hbsize=1e4,dt=0.01)
times=seq(0,365*2,1)

tempFunc=function(time,deltaT){
  T=15*(1-cos(2*pi*(time+80)/365))+deltaT
  return(T)}

gFunc=function(T){
  v=(1/60)*(1-((T-20)/20)^2)
  return(max(0,v))}

ccFunctions<-list(
                         reproFunc=function(x,time,species,strain){
                           A0=600;  q=11.5
                           reprod=q*x$Nematodes['adults',1]*exp(-x$Nematodes['adults',1]/A0)
                           return(max(0,reprod))
                         },
                         deathFunc=function(stage,x,time,species,strain){
                           a=c(0.05,0.05)
                           v=a[stage]
                           return(max(0,v))
                         },
                         develFunc=function(stage,x,time,species,strain){
                           T=tempFunc(time,deltaT)
                           v=gFunc(T)
                           return(v)
                         },
                         durationFunc=function(stage,x,time,species,strain){
                           if (time==0){
                             T=tempFunc(time,deltaT)
                             v=1/ccFunctions$develFunc(stage,x,time,species,strain)}
                           return(v)
                         },
                         immigrationFunc=function(stage,x,time,species,strain){
                           v=0
                           if (stage==2){if (time>=0 & time <=0.1){v=1}}
                           return(v)
                         },
                         emigrationFunc=function(stage,x,time,species,strain){return(0)}
)


Adults=matrix(0,nrow=length(times),ncol=2)
Juvs=matrix(0,nrow=length(times),ncol=2)
Taus=matrix(0,nrow=length(times),ncol=2)
ct=0
for (deltaT in c(0,3)){

  ct=ct+1

  modelOutput=popModel(
    numSpecies=1,
    numStages=2,
    timeDependLoss=FALSE,
    timeDependDuration=TRUE,
    ICs=list(matrix(0,nrow=2,ncol=1)),
    timeVec=times,
    solverOptions=solver.options,
    rateFunctions=ccFunctions,
    stageNames=list(c('juveniles','adults')),
    speciesNames=c('Nematodes'),
    saveFig=TRUE,
    figType='png'
    )

  Adults[,ct]=modelOutput[,'Nematodes.adults']
  Juvs[,ct]=modelOutput[,'Nematodes.juveniles']
  Taus[,ct]=modelOutput[,'dur.Nematodes.juveniles']
  
}

dev.new(bg='white')
par(mar=c(5,5,2,2))
T=seq(0,40,0.5); g=T*0
for (i in seq(1,length(T))){g[i]=gFunc(T[i])}
plot(T,g,type='l',xlab='Temperature (oC)',ylab='Juvenile development rate (/d)',cex.axis=1.5,cex.lab=1.5)
dev.copy2eps(file='EnvTempDev.eps')
#dev.print(png,'EnvTempDev.png',res=100,width=5,height=5, units='in')

dev.new(bg='white')
par(mar=c(5,5,2,2))
time=modelOutput[,'time']
temps=tempFunc(time,0)
plot(c(0,max(times)),c(0,35),type="n",xlab='Time',ylab='Temperature',cex.axis=1.75,cex.lab=2)
lines(time,temps,col=1,lty=1,lwd=2)
temps=tempFunc(time,3)
lines(time,temps,col=2,lty=1,lwd=2)
legend('top',c('dT=0','dT=3'),lty=1,col=seq(1,2),lwd=2,cex=1.25,bty='n')
dev.copy2eps(file='EnvTemp.eps')
#dev.print(png,'EnvTemp.png',res=100,width=5,height=5, units='in')


dev.new(bg='white')
par(mar=c(5,5,2,2))
plot(c(0,2*365),c(60,200),type="n",xlab='Time',ylab='Juvenile Stage Duration (d)',cex.axis=1.75,cex.lab=2)
lines(time,Taus[,1],col=1,lty=1,lwd=2)
lines(time,Taus[,2],col=2,lty=1,lwd=2)
legend('top',c('dT=0','dT=3'),lty=1,col=seq(1,2),lwd=2,cex=1.25,bty='n')
dev.copy2eps(file='EnvTau.eps')
#dev.print(png,'EnvTau.png',res=100,width=5,height=5, units='in')

dev.new(bg='white')
par(mar=c(5,5,2,2))
plot(c(0,2*365),c(0,100),type="n",xlab='Time',ylab='Adult Density',cex.axis=1.75,cex.lab=2)
lines(time,Adults[,1],col=1,lty=1,lwd=2)
lines(time,Adults[,2],col=2,lty=1,lwd=2)
legend('top',c('dT=0','dT=3'),lty=1,col=seq(1,2),lwd=2,cex=1.25,bty='n')
dev.copy2eps(file='EnvAdults.eps')
#dev.print(png,'EnvTau.png',res=100,width=5,height=5, units='in')




 

