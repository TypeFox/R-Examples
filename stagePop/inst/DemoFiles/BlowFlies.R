#Nicholson's Blow Flies (Gurney et al 1983)
#Reproducing the first example in 'The systematic formulation of tractable single-species population models incorporating age structure'
#by WSC Gurney, RM Nisbet and J Lawton, Journal of Animal Ecology, 52, p479-495, 1983

#library(stagePop)

#All the vectors are specified in the order of the life cycle
#e.g. start with eggs and finish with reproducing adults

#solver.options=list(DDEsolver='deSolve',atol=1e-3,rtol=1e-3,hbsize=1e4)
solver.options=list(DDEsolver='PBS',tol=1e-4,hbsize=1e4,dt=0.01)

#--------------------------DEFINE RATE FUNCTIONS---------------------------------------
blowFliesFunctions <- list(
                           reproFunc=function(x,time,species,strain){
                             A0=600;  q=8.5
                             reprod=q*x$blowflies['adults',1]*exp(-x$blowflies['adults',1]/A0)
                             return(reprod)
                           },
                           deathFunc=function(stage,x,time,species,strain){
                           #per capita death rate (/d)
                             a=c(0.07,0.004,0.003,0.0025,0.27)
                             return(a[stage])
                           },
                           durationFunc=function(stage,x,time,species,strain){
                           #duration of each stage in days
                             a=c(0.6,5.0,5.9,4.1)
                             return(a[stage])
                           },
                           immigrationFunc=function(stage,x,time,species,strain){
                             v=0
                             if (stage==5){if (time>=0 & time<=1){v=100}}
                             return(v)
                           },
                           emigrationFunc=function(stage,x,time,species,strain){return(0)}
                           )

#-------------------------------RUN MODEL---------------------------------------

modelOutput = popModel(
  numSpecies=1,
  numStages=5,
  ICs=list(matrix(0,nrow=5,ncol=1)),
  timeVec=seq(0,200,0.5),
  timeDependLoss=TRUE,
  timeDependDuration=FALSE,
  rateFunctions=blowFliesFunctions,
  solverOptions=solver.options,
  stageNames=list(c('eggs','larvae','pupae','juveniles','adults')),
  speciesNames=c('blowflies'),
  saveFig=TRUE,
  figType='png'
  )





