initConditions=function(definePop,x,rateFunctions,startTime,speciesNames,stageNames){

  #convert the input ICs to the vector needed by the derivs function
  #x is named ICs list
  
#  print('enter initConditions')

  allInitN={};  allInitP={};  allInitTau={}
  initP={}; initTau={}
  numSpecies=definePop$numSpecies

  for (species in seq(1,numSpecies)){

    initN=x[[species]]

    if (is.matrix(initN)){initN=as.vector(initN)} 

    timeDependLoss=definePop$timeDependLoss[species]
    timeDependDuration=definePop$timeDependDuration[species]
    numStages=definePop$numStages[species]
    numStrains=definePop$numStrains[species]

    if (numStages>1){

      Tau0=seq(1,(numStages-1)*numStrains)
      ct=1
      for (r in seq(1,numStrains)){
        for (i in seq(1,numStages-1)){
          Tau0[ct]=rateFunctions$durationFunc(i,x,startTime,species,r)
          ct=ct+1
        }
      }

      
      if ((timeDependLoss | timeDependDuration) & numStages>1){
        initP=0*seq(1,(numStages-1)*numStrains)
        ct=1
        for (r in seq(1,numStrains)){
          for (i in seq(1,numStages-1)){
            initP[ct]=exp(-rateFunctions$deathFunc(i,x,startTime,species,r)*rateFunctions$durationFunc(i,x,startTime,species,r))
            ct=ct+1
          }
        }
      }else{initP={}}
      
      if (timeDependDuration){
        initTau=Tau0
      }else{initTau={}}
    }else{initP={}; initTau={}}
    
    allInitN=c(allInitN,initN)
    allInitP=c(allInitP,initP)
    allInitTau=c(allInitTau,initTau)
  }


  return(c(allInitN,allInitP,allInitTau))
}
