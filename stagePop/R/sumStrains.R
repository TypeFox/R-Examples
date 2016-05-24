#' sumStrains
#' @param namedOutput The output martix from popModel
#' @param numSpecies Number of species (integer)
#' @param numStages Number of life-stages in each species (vector)
#' @param numStrains Number of strains in each species (vector)
#' @param speciesNames A vector of strings containing the name for each species. This is only used for plotting purposes. 
#' @param stageNames A list of n vectors (where n is the number of species) containing the names of each stage for each species. 
#' @param timeDependLoss A vector specifying TRUE/FALSE for each species. If a species has any time dependent per capita death rates (e.g. density dependent death rates) this is TRUE. 
#' @param timeDependDuration A vector specifying TRUE/FALSE for each species. If a species has any time dependent stage durations it is TRUE for that species. 
#' @export

sumStrains=function(namedOutput,numSpecies,numStages,numStrains,speciesNames,stageNames,timeDependLoss,timeDependDuration){


  #used the fully named modelOutput matrix and then sum over the rows using the names and grepl
  newOut=NA*namedOutput
  newNames=rep('NA',length(namedOutput[1,]))

  #time column first
  newOut[,1]=namedOutput[,1]
  newNames[1]='time'

  #state variables
  ct=2  
  for (i in seq(1,numSpecies)){
    for (j in seq(1,numStages[i])){
      text=paste(speciesNames[i],'.',stageNames[[i]][j],sep='')
      key=paste('^',text,sep='')
      mat=namedOutput[,grepl(key,colnames(namedOutput))]
      if (is.matrix(mat)){newOut[,ct]=rowSums(mat)}else{newOut[,ct]=mat}
      newNames[ct]=text
      ct=ct+1
    }}
  #survival probabilities
  for (i in seq(1,numSpecies)){
    for (j in seq(1,(numStages[i]-1))){
      if ((timeDependLoss[i] | timeDependDuration[i]) & numStages[i]>1){
        text=paste('prob.',speciesNames[i],'.',stageNames[[i]][j],sep='')
        key=paste('^',text,sep='')
        mat=namedOutput[,grepl(key,colnames(namedOutput))]
        if (is.matrix(mat)){newOut[,ct]=rowMeans(mat)}else{newOut[,ct]=mat}
        newNames[ct]=text
        ct=ct+1
      }
    }}
 #stage durations
  for (i in seq(1,numSpecies)){
    for (j in seq(1,(numStages[i]-1))){
      if (timeDependDuration[i] & numStages[i]>1){
        text=paste('dur.',speciesNames[i],'.',stageNames[[i]][j],sep='')
        key=paste('^',text,sep='')
        mat=namedOutput[,grepl(key,colnames(namedOutput))]
        if (is.matrix(mat)){newOut[,ct]=rowMeans(mat)}else{newOut[,ct]=mat}
        newNames[ct]=text
        ct=ct+1
      }
    }}
  #rates
  for (i in seq(1,numSpecies)){
    for (j in seq(1,(numStages[i]))){
      text=paste('dot.',speciesNames[i],'.',stageNames[[i]][j],sep='')
      key=paste('^',text,sep='')
      mat=namedOutput[,grepl(key,colnames(namedOutput))]
      if (is.matrix(mat)){newOut[,ct]=rowSums(mat)}else{newOut[,ct]=mat}
      newNames[ct]=text
      ct=ct+1
    }}

  
  newOut=newOut[,1:(ct-1)]
  colnames(newOut)=newNames[1:(ct-1)]
  
  return(newOut)  
}
