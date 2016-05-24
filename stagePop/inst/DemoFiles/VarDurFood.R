#Variable stage duration (Nisbet and Gurney)
#used scaled version as in N&G appendix 2

#Food (F) is eaten by larvae (L)
#when larvae reach a certain mass they become adults (A).
#variables are scaled
#first 'species' is food

library(stagePop)

#All the vectors are specified in the order of the life cycle
#e.g. start with eggs and finish with reproducing adults

#PARAMETERS
fs=1     #rate of food supply
fmax=3   #max rate of food supply
F0=0.1   #initial food density
m=1      #number of mass units a larva must increase by before becoming an adult
epsilon=1#constant of proportionality between development and food consumption
K=1      #half sat constant for food consumption
q=5#500      #reproduction rate
dA=2#10     #adult death rate
dL=log(q/dA)  #larval death rate

solver.options=list(DDEsolver='deSolve',atol=1e-9,rtol=1e-9,method='lsoda',hbsize=1e4)
#solver.options=list(DDEsolver='PBS',tol=1e-9,hbsize=1e4,dt=0.01)

#--------------------------DEFINE RATE FUNCTIONS---------------------------------------
varDurFoodFunctions <- list(
                            reproFunc=function(x,time,species,strain){
                              if (species==1){reprod=fs}
                              if (species==2){reprod=q*x$Damselfly['adults',1]}
                              return(max(0,reprod))
                            },
                            deathFunc=function(stage,x,time,species,strain){
                              if(species==1){v=fmax*x$Damselfly['larvae',1]/(K+x$Food[1,1])}
                              if(species==2){a=c(dL,dA);v=a[stage]}
                              return(max(0,v))
                            },
                            durationFunc=function(stage,x,time,species,strain){
                              if (time==0 & species==2 & stage==1){v=m/(epsilon*fmax*F0/(K+F0))}
                              return(v)
                            },
                            develFunc=function(stage,x,time,species,strain){
                              if (species==2 & stage==1){v=epsilon*fmax*x$Food[1,1]/(K+x$Food[1,1])}
                              return(v)
                            },
                            immigrationFunc=function(stage,x,time,species,strain){
                              v=0
                              if (species==2 & stage==1){if (time>=0 & time <=0.1){v=1}}
                              return(v)
                            },
                            emigrationFunc=function(stage,x,time,species,strain){return(0)}
                            )


#-------------------------------RUN MODEL---------------------------------------
  
modelOutput=popModel(
  numSpecies=2,
  speciesNames=c('Food','Damselfly'),
  numStages=c(1,2),
  stageNames=list('one',c('larvae','adults')),
  numStrains=c(1,1),
  timeDependLoss=c(TRUE,FALSE),
  timeDependDuration=c(FALSE,TRUE),
  ICs=list(matrix(F0,1,1),matrix(0,nrow=2,ncol=1)),
  timeVec=seq(0,30,0.1),
  solverOptions=solver.options,
  rateFunctions=varDurFoodFunctions,
  sumOverStrains=TRUE
)

tauStar=(1/dL)*log(q/dA)
qdash=q*tauStar
dAdash=dA*tauStar
fmaxdash=epsilon*fmax*tauStar/m
Kdash=K/(fs*tauStar)
print(c(qdash,dAdash,fmaxdash,Kdash))

tauStar=log(q/dA)/dL
Lstar=fs*tauStar*epsilon/m
Fstar=K/((epsilon*fmax*tauStar/m)-1)
Astar=dL*Lstar/(q-dA)

print(c('tauStar','Fstar','Lstar','Astar'))
print(c(tauStar,Fstar,Lstar,Astar))
print(modelOutput[length(modelOutput[,1]),c('dur.Damselfly.larvae','Food.one','Damselfly.larvae','Damselfly.adults')])

#plot larval stage duration
dev.new()
par(mar=c(5,5,1,1))
plot(modelOutput[,'time'],modelOutput[,'dur.Damselfly.larvae'],type='l',xlab='time',ylab='Larval stage duration',cex.axis=1.5,cex.lab=2.5,lwd=2)
#dev.copy2eps(file='Foodtau.eps')

#plot rate of uptake of food 
foodUptake=epsilon*fmax*modelOutput[,'Food.one']/(K+modelOutput[,'Food.one'])
dev.new()
par(mar=c(5,5,1,1))
plot(modelOutput[,'time'],foodUptake,type='l',xlab='time',ylab='Food uptake rate',cex.axis=1.5,cex.lab=2.5,lwd=2)
#dev.copy2eps(file='FoodUptake.eps')

