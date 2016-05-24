addNamesToList=function(matList,speciesNames,stageNames,numSpecies){

  newList=list()
  for (i in seq(1,numSpecies)){
    mat=as.matrix(matList[[i]])
    rownames(mat)=stageNames[[i]]
    newList=c(newList,list(mat))
  }
  names(newList)=speciesNames
  return(newList)
  
}
