namingVariables=function(numSpecies,numStages,numStrains,speciesNames,stageNames,timeDependLoss,timeDependDuration){

#  print('enter namingVariables')
  #produces a vector containing the names of the columns in the output from the DDE solvers

  L=sum(numStages*numStrains)
  stVarNames=seq(1,L)
  probNames=NA*seq(1,L)
  durNames=NA*seq(1,L)
  dotNames=seq(1,L)

  
  ct=1
  for (i in seq(1,numSpecies)){
    species.name=speciesNames[i]
    for (k in seq(1,numStrains[i])){
      if (numStrains[i]==1){strain.name=NULL}else{strain.name=paste('.strain',k,sep='')}
      for (j in seq(1,numStages[i])){
        stage.name=stageNames[[i]][j]
        stVarNames[ct]=paste(species.name,'.',stage.name,strain.name,sep='')
        ct=ct+1
      }}}
  
  ct=1
    for (i in seq(1,numSpecies)){
      species.name=speciesNames[i]
      if ((timeDependLoss[i] | timeDependDuration[i]) & numStages[i]>1){
        for (k in seq(1,numStrains[i])){
          if (numStrains[i]==1){strain.name=NULL}else{strain.name=paste('.strain',k,sep='')}
          for (j in seq(1,(numStages[i]-1))){
            stage.name=stageNames[[i]][j]
            probNames[ct]=paste('prob.',species.name,'.',stage.name,strain.name,sep='')
            ct=ct+1
          }}}}

  ct=1
  for (i in seq(1,numSpecies)){
    species.name=speciesNames[i]
    if (timeDependDuration[i] & numStages[i]>1){
        for (k in seq(1,numStrains[i])){
          if (numStrains[i]==1){strain.name=NULL}else{strain.name=paste('.strain',k,sep='')}
          for (j in seq(1,(numStages[i]-1))){
            stage.name=stageNames[[i]][j]
            durNames[ct]=paste('dur.',species.name,'.',stage.name,strain.name,sep='')
            ct=ct+1
          }}}}

  ct=1
  for (i in seq(1,numSpecies)){
    species.name=speciesNames[i]
    for (k in seq(1,numStrains[i])){
      if (numStrains[i]==1){strain.name=NULL}else{strain.name=paste('.strain',k,sep='')}
      for (j in seq(1,numStages[i])){
        stage.name=stageNames[[i]][j]
        dotNames[ct]=paste('dot.',species.name,'.',stage.name,strain.name,sep='')
        ct=ct+1
      }}}

  vNames=c(stVarNames,probNames[!is.na(probNames)],durNames[!is.na(durNames)],dotNames)

  return(vNames)
  }
