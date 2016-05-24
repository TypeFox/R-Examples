#This classic predator prey model (with no stage structure)

library(stagePop)

growthRatePred=10
growthRatePrey=1
deathRatePrey=1
deathRatePred=0.5


#solver.options=list(DDEsolver='deSolve',atol=1e-6,rtol=1e-6,method='lsoda',hbsize=1e4)
solver.options=list(DDEsolver='PBS',tol=1e-7,hbsize=1e4,dt=0.01)

ppFunctions <- list(
                    reproFunc=function(x,time,species,strain){
                      if (species==1){reprod=growthRatePrey*x$prey['adults',1]}
                      if(species==2){reprod=growthRatePred*x$prey['adults',1]*x$predator['adults',1]}
                      return(max(0,reprod))
                    },
                    deathFunc=function(stage,x,time,species,strain){
                      if (species==1){v=deathRatePrey*x$predator['adults',1]}
                      if (species==2){v=deathRatePred}
                      return(max(0,v))
                    },
                    durationFunc=function(stage,x,time,species,strain){
                      return(0)
                    },
                    immigrationFunc=function(stage,x,time,species,strain){
                      return(0)
                    },
                    emigrationFunc=function(stage,x,time,species,strain){
                      return(0)
                    }
                    )
#-------------------------------RUN MODEL----------------------------------------------------
modelOutput=popModel(
  numSpecies=2,
  numStages=c(1,1),
  timeDependLoss=c(TRUE,FALSE),
  timeDependDuration=c(FALSE,FALSE),
  ICs=list(matrix(0.3,1,1),matrix(1,1,1)),
  timeVec=seq(0,100,0.1),
  solverOptions=solver.options,
  plotFigs=TRUE,
  rateFunctions=ppFunctions,
  speciesNames=c('prey','predator'),
  stageNames=list('adults','adults')
)

dev.new()
par(mar=c(5,5,2,2))
plot(modelOutput[,2],modelOutput[,3],type='l',xlab='prey',ylab='predator',cex.lab=2,cex.axis=1.5)
#dev.copy2eps(file='predprey1C.eps')

x=0.3; y=1
C=growthRatePrey*log(y)-deathRatePrey*y-growthRatePred*x+deathRatePred*log(x)

x=modelOutput[,2]; y=modelOutput[,3]
Cmod=growthRatePrey*log(y)-deathRatePrey*y-growthRatePred*x+deathRatePred*log(x)
print(c(C,summary(Cmod)))





