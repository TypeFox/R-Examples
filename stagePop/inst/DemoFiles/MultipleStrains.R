library(stagePop)

#solver.options=list(DDEsolver='deSolve',atol=1e-5,rtol=1e-7,method='lsoda',hbsize=1e4)
solver.options=list(DDEsolver='PBS',tol=1e-7,hbsize=1e4,dt=0.01)

#parameters
Rin=10; V=1; K=1
num.strains=6
Yield=0.5#g of bacteria produced for 1 g of food
if (num.strains>1){Gmax=2+seq(1,num.strains)}else{Gmax=2}

case=2
if (case==1){
  num.stages=1;  stage.names='reproductive';  start=0.1;numDays=30}
if (case==2){
  num.stages=2;  stage.names=c('lagged','reproductive');  start=c(0,0.1);numDays=100}
  

strainsFunctions <- list(
                           reproFunc=function(x,time,species,strain){
                             if (species==1){reprod=Rin*V}
                             if (species==2){
                               reprod=x$Bacteria['reproductive',strain]*
                                 Gmax[strain]*x$Resource['food',1]/(x$Resource['food',1]+K)}
                             return(reprod)
                           },
                           deathFunc=function(stage,x,time,species,strain){
                           #per capita death rate (/d)
                             if (species==1){
                               uptake=0*seq(1,num.strains)
                               for (s in seq(1,num.strains)){
                                 uptake[s]=(Gmax[s]/(x$Resource['food',1]+K))*(x$Bacteria['reproductive',s]/Yield)
                               }
                               death=sum(uptake)+V
                             }
                             if (species==2){
                               if (stage==1){if(num.stages==2){death=0}else{death=V}}
                               if (stage==2){death=V}
                             }
                             return(death)
                           },
                           durationFunc=function(stage,x,time,species,strain){
                           #duration of each stage in days
                             durations=2*seq(1,num.strains)
                             return(durations[strain])
                           },
                           immigrationFunc=function(stage,x,time,species,strain){return(0)},
                           emigrationFunc=function(stage,x,time,species,strain){return(0)}
                           )


modelOutput = popModel(
  numSpecies=2,
  numStrains=c(1,num.strains),
  numStages=c(1,num.stages),
  ICs=list(matrix(Rin,nrow=1,ncol=1),matrix(start,nrow=num.stages,ncol=num.strains)),
  timeVec=seq(0,numDays,0.5),
  timeDependLoss=c(TRUE,FALSE),
  timeDependDuration=c(FALSE,FALSE),
  rateFunctions=strainsFunctions,
  solverOptions=solver.options,
  stageNames=list(c('food'),stage.names),
  speciesNames=c('Resource','Bacteria'),
  sumOverStrains=FALSE,
  plotStrainsFig=TRUE
  )


