vecToMatList=function(vec,numSpecies,numStages,numStrains,speciesNames,stageNames){
  #convert the vector used by the deriv functions to a list of matrices for use in the rate functions
  #each matrix has row names which are the stage names for that species
  #each matrix has a species name
  start=1
  matList=list()
  for (i in seq(1,numSpecies)){
    len=numStages[i]*numStrains[i]
    speciesVec=vec[start:(start+len-1)]
    speciesMat=matrix(speciesVec,nrow=numStages[i],ncol=numStrains[i])
    rownames(speciesMat)=stageNames[[i]]
    matList=c(matList,list(speciesMat))
    start=start+len
  }
  names(matList)=speciesNames
  return(matList)
}



