#Reproducing the second example in 'The systematic formulation of tractable single-species population models incorporating age structure'
#by WSC Gurney, RM Nisbet and J Lawton, Journal of Animal Ecology, 52, p479-495. 1983

library(stagePop)


#All the vectors are specified in the order of the life cycle
#e.g. start with eggs and finish with reproducing adults

#solver.options=list(DDEsolver='deSolve',atol=1e-4,rtol=1e-4,method='lsoda',hbsize=1e6)
solver.options=list(DDEsolver='PBS',tol=1e-6,hbsize=1e3,dt=0.1)

case=2#choose case (1 or 2)
if (case==1){
  num.stages=2
  stage.names=c('larvae','adults')
}else{
  num.stages=3
  stage.names=c('larvae','adults','dead adults')}

#--------------------------DEFINE RATE FUNCTIONS---------------------------------------
larvalCompFunctions <- list(
                            reproFunc=function(x,time,species,strain){
                              reprod=9.4*x$flies['adults',1]
                              return(reprod)
                            },
#-----------------------------------------------------------
                            deathFunc=function(stage,x,time,species,strain){
                              if (stage==1){v=5e-5*x$flies['larvae',1]}
                              if (stage>=2){if (case==1){v=0.2}else{v=0}}
                              return(v)
                            },
#-----------------------------------------------------------
                            durationFunc=function(stage,x,time,species,strain){
                              a=c(28,5)
                              return(a[stage])
                            },
#-----------------------------------------------------------
                            immigrationFunc=function(stage,x,time,species,strain){
                              v=0
                              if (time>=0 & time<=1){
                                if (stage==2){v=20}}
                              return(v)
                            },
                            emigrationFunc=function(stage,x,time,species,strain){return(0)}
                            )

#-------------------------------RUN MODEL---------------------------------------

modelOutput=popModel(
  numSpecies=1,
  numStages=num.stages,
  ICs=list(matrix(0,nrow=num.stages,ncol=1)),
  timeVec=seq(0,500,0.5),
  timeDependLoss=TRUE,
  solverOptions=solver.options,
  rateFunctions=larvalCompFunctions,
  stageNames=list(stage.names),
  speciesNames='flies'
  )







