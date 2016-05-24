#Example where the juvenile stage duration is dependent upon temperature.

library(stagePop)

#solver.options=list(DDEsolver='deSolve',atol=1e-6,rtol=1e-6,hbsize=1e5)
solver.options=list(DDEsolver='PBS',tol=1e-8,hbsize=1e4,dt=0.01)

maxDur=200
minDur=60
tempFunc=function(time){
  T=15*(1-cos(2*pi*(time+80)/365))
  return(T)}

tauFunc=function(T){
  v=min(minDur+((T-20)/2)^2,maxDur)
  return(v)}

varDurEnvFunctions<-list(
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
                           T=tempFunc(time)
                           v=1/tauFunc(T)
                           return(v)
                         },
                         durationFunc=function(stage,x,time,species,strain){
                           if (time==0){
                             T=tempFunc(time)
                             v=tauFunc(T)}
                           return(v)
                         },
                         immigrationFunc=function(stage,x,time,species,strain){
                           v=0
                           if (stage==2){if (time>=0 & time <=0.1){v=1}}
                           return(v)
                         },
                         emigrationFunc=function(stage,x,time,species,strain){return(0)}
)


modelOutput=popModel(
  numSpecies=1,
  numStages=2,
  timeDependLoss=FALSE,
  timeDependDuration=TRUE,
  ICs=list(matrix(0,nrow=2,ncol=1)),
  timeVec=seq(0,365*6,1),
  solverOptions=solver.options,
  rateFunctions=varDurEnvFunctions,
  stageNames=list(c('juveniles','adults')),
  speciesNames=c('Nematodes')
  )



dev.new()
par(mar=c(5,5,5,5))
T=seq(0,40,0.5); Tau=T*0
for (i in seq(1,length(T))){Tau[i]=tauFunc(T[i])}
plot(T,Tau,type='l',xlab='Temperature (oC)',ylab='Juvenile stage Duration (days)',cex.axis=1.75,cex.lab=2)
#dev.copy2eps(file='EnvTempTau.eps')

dev.new()
par(mar=c(5,5,5,5))
time=modelOutput[,'time']
temps=tempFunc(time)
plot(c(0,2*365),c(0,max(temps)),type="n",xlab='Time',ylab='Temperature',cex.axis=1.75,cex.lab=2)
lines(time,temps,col=1,lty=1,lwd=2)
#dev.copy2eps(file='EnvTemp.eps')


dev.new()
par(mar=c(5,5,5,5))
tau=modelOutput[,'dur.Nematodes.juveniles']
tempsTau=0*time
for (i in seq(1,length(time))){tempsTau[i]=tauFunc(tempFunc(time[i]))}
plot(c(0,2*365),c(minDur*0.9,0.9*maxDur),type="n",xlab='Time',ylab='Juvenile Stage Duration (d)',cex.axis=1.75,cex.lab=2)
lines(time,tau,col=1,lty=1,lwd=2)
lines(time,tempsTau,col=2,lty=1,lwd=2)
legend('topleft',c('from StagePop','if T(t)=Tc'),lty=1,col=seq(1,2),lwd=2,cex=1.75)
#dev.copy2eps(file='EnvTau.eps')


