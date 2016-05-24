#This classic predator prey model but with stage structure for the predator
#only adult predators eat prey

library(stagePop)

growthRatePred=10
growthRatePrey=1
deathRatePrey=1
deathRatePred=c(1.0,0.5)
preyCarryCapacity=1
  
#solver.options=list(DDEsolver='deSolve',atol=1e-6,rtol=1e-6,method='lsoda',hbsize=1e4)
solver.options=list(DDEsolver='PBS',tol=1e-7,hbsize=1e4,dt=0.01)

case=1; print(paste('case',case))
if (case==1){juvPredDuration=0.1;lenTime=100}
if (case==2){juvPredDuration=0.1;lenTime=300}
if (case==3){juvPredDuration=1.8;lenTime=300}
if (case==4){juvPredDuration=0.1;lenTime=400;deathRatePred[1]=0}
if (case==5){juvPredDuration=5;lenTime=400;deathRatePred[1]=0}
if (case==6){juvPredDuration=15;lenTime=400;deathRatePred[1]=0}
if (case==7){juvPredDuration=20;lenTime=400;deathRatePred[1]=0}

ppFunctions <- list(
                    reproFunc=function(x,time,species,strain){
                      if (species==1){reprod=growthRatePrey*x$Prey['adult',1]}
                      if (species==2){reprod=growthRatePred*x$Prey['adult',1]*x$Predator['adult',1]}
                      return(max(0,reprod))
                    },
                    deathFunc=function(stage,x,time,species,strain){
                      if (species==1){
                        if (case==1){v=deathRatePrey*x$Predator['adult',1]}
                        if (case>1){v=deathRatePrey*x$Predator['adult',1]+growthRatePrey*x$Prey['adult',1]/preyCarryCapacity}}
                      if (species==2){v=deathRatePred[stage]}
                      return(max(0,v))
                    },
                    durationFunc=function(stage,x,time,species,strain){
                      return(juvPredDuration)
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
  numStages=c(1,2),
  timeDependLoss=c(TRUE,FALSE),
  timeDependDuration=c(FALSE,FALSE),
  ICs=list(matrix(0.3),matrix(c(0,1),nrow=2,ncol=1)),
  timeVec=seq(0,lenTime,0.1),
  solverOptions=solver.options,
  plotFigs=TRUE,
  rateFunctions=ppFunctions,
  stageNames=list('adult',c('juvenile','adult')),
  speciesNames=c('Prey','Predator')
  )






