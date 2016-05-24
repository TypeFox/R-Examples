#Host parasitoid model 
#Briggs CJ. 1993. Competition among parasitoid species on a stage-structured host and its effect on host suppression. The American Naturalist 141, 372-397.

library(stagePop)

#All the vectors are specified in the order of the life cycle
#e.g. start with eggs and finish with reproducing adults
attackRateP=1
attackRateQ=2
TE=0.5
TL=0.5
TJP=0.4
TJQ=0.4
deathE=0.1; deathL=0.1; deathA=0.1
deathJP=0.1; deathP=8.0
deathJQ=0.1; deathQ=8.0
rho=33 #total lifetime fecundity 

#Analytical Solution#---------------------------------------------
BriggAnalytic=function(Ppresent,attackRates,durations,deathRates,rho){

  attackRateP=attackRates[1];  attackRateQ=attackRates[2]
  TE=durations[1];  TL=durations[2]
  TJP=durations[3];  TJQ=durations[4]
  deathE=deathRates[1]; deathL=deathRates[2]; deathA=deathRates[3]
  deathJP=deathRates[4]; deathP=deathRates[5]
  deathJQ=deathRates[6]; deathQ=deathRates[7]

  if (Ppresent){
    Pstar=(log(rho)-deathL*TL-deathE*TE)/(attackRateP*TE)
    Qstar=0
  }else{
    Qstar=(log(rho)-deathL*TL-deathE*TE)/(attackRateQ*TL)
    Pstar=0
  }
  
  yP=attackRateP*Pstar+deathE
  yQ=attackRateQ*Qstar+deathL
  
  if (Ppresent){
    Estar=deathP/(attackRateP*exp(-deathJP*TJP))
    Astar=Estar*yP/(rho*deathA*(1-exp(-yP*TE)))
    Lstar=(rho*deathA*Astar*(1-exp(-yQ*TL))*exp(-yP*TE))/yQ
  }else{
    Lstar=deathQ/(attackRateQ*exp(-deathJQ*TJQ))
    Astar=Lstar*yQ/(rho*deathA*(1-exp(-yQ*TL))*exp(-yP*TE))
    Estar=Astar*rho*deathA*(1-exp(-yP*TE))/yP
  }
  
  return(c(Estar,Lstar,Astar,Pstar,Qstar))
}

#-------------------------------------------------------------
durations=c(c(TE,TL),TJP,TJQ)
deathRates=c(deathE,deathL,deathA,deathJP,deathP,deathJQ,deathQ)
attackRates=c(attackRateP,attackRateQ)
analy=BriggAnalytic(FALSE,attackRates,durations,deathRates,rho)
LstarQ=analy[2]
AstarQ=analy[3]
Qstar=analy[5]

BriggsFunctions <- list(
              reproFunc=function(x,time,species,strain){
                if (species==1){reprod=rho*deathA*x$Host['adults',1]}
                if (species==2){reprod=attackRateP*x$'Egg Parasitoid'['adults',1]*x$Host['eggs',1]}
                if (species==3){reprod=attackRateQ*x$'Larval Parasitoid'['adults',1]*x$Host['larvae',1]}
                return(max(0,reprod))
              },
              deathFunc=function(stage,x,time,species,strain){
                if (species==1){a=c(deathE,deathL,deathA);v=a[stage]
                    if (stage==1){v=a[stage]+attackRateP*max(x$'Egg Parasitoid'['adults',1],0)}
                    if (stage==2){v=a[stage]+attackRateQ*max(x$'Larval Parasitoid'['adults',1],0)}}
                if (species==2){a=c(deathJP,deathP);v=a[stage]}
                if (species==3){a=c(deathJQ,deathQ);v=a[stage]}
                return(max(0,v))
                },
               durationFunc=function(stage,x,time,species,strain){
                 if (species==1){a=c(TE,TL)}
                 if (species==2){a=TJP}
                 if (species==3){a=TJQ}
                 return(a[stage])
               },
               immigrationFunc=function(stage,x,time,species,strain){
                 v=0
                 if (species==1){if (time>=0 & time<=0.1){f1=rho*deathA*AstarQ
                   if (stage==1){v=f1}
                   if (stage==2){v=f1*exp(-deathE*TE)}
                   if (stage==3){v=f1*exp(-deathE*TE-deathL*TL)}}}
                 if(species==2){if(time>=20 & time<=20.1){if(stage==2){v=1}}}
                 if(species==3){if(time>=0 & time<=0.1){if(stage==1){v=attackRateQ*Qstar*LstarQ}}}
                 return(v)
               },
               emigrationFunc=function(stage,x,time,species,strain){return(0)}
               )

modelOutput=popModel(
  numSpecies=3,
  numStages=c(3,2,2),
  timeVec=seq(0,50,0.1),
  rateFunctions=BriggsFunctions,
  timeDependLoss=c(TRUE,FALSE,FALSE),
  timeDependDuration=c(FALSE,FALSE,FALSE),
  ICs=list(matrix(0,nrow=3,ncol=1),matrix(0,nrow=2,ncol=1),matrix(0,nrow=2,ncol=1)),
  solverOptions=list(DDEsolver='PBS',tol=1e-7,hbsize=1e4,dt=0.01),
  speciesNames=c('Host','Egg Parasitoid','Larval Parasitoid'),
  stageNames=list(c('eggs','larvae','adults'),c('eggs','adults'),c('eggs','adults')),
  saveFig=TRUE,
  figType='png'
  )


#check results against equilibrium values (Pg 392 in Briggs 1993)
time=modelOutput[,1]
print('model')
L=length(modelOutput[,1])
print(round(modelOutput[L-1,c(2:4,6,8)],3))
if (modelOutput[L-1,8]>modelOutput[L-1,6]){
  print('Just Q')
  Ppresent=FALSE
}else{
  print('Just P')
  Ppresent=TRUE  }

print('Analytical')
print(round(BriggAnalytic(Ppresent,attackRates,durations,deathRates,rho),3))



