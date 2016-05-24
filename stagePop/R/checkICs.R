checkICs=function(ICs,numSpecies,numStages,numStrains){

#  print('enter checkICs')
  
  if (!is.list(ICs)){
    stop("Error in ICs: The ICs must be a list containing matrices for each species. For species i the matrix must have nrow=numStages[i] and ncol=numStrains[i].")
  }else{
    for (i in seq(1,numSpecies)){
      if (!is.matrix(ICs[[i]])){
        stop(paste('Error in ICs: The entry in ICs for species',i,'is not a matrix. The ICs must be a list containing matrices for each species. For species i the matrix must have nrow=numStages[i] and ncol=numStrains[i]. Use the matrix function e.g. matrix(0,nrow=1,ncol=1) function.'))
      }else{
        if (ncol(ICs[[i]])!=numStrains[i]){stop(paste('Error in ICs: The number of columns in the matrix for species',i,'in the ICs must be equal to the number of strains specified for species',i,'in numStrains. In this case ncol=',ncol(ICs[[i]]),'and numStrains is',numStrains[i]))}
        if (nrow(ICs[[i]])!=numStages[i]){stop(paste('Error in ICs: The number of rows in the matrix for species',i,'in the ICs must be equal to the number of stages specified for species',i,'in numStages. In this case nrow=',nrow(ICs[[i]]),'and numStages is',numStages[i]))}
        #check zeros for all but last (and only) stage
        if (numStages[i]>1){
          for (k in seq(1,numStrains[i])){
            if (sum(ICs[[i]][1:(numStages[i]-1),k])>0) stop(paste('Error in ICs: Non-zero initial conditions are only allowed for single stage species or for the last stage of a multi-stage species - please change ICs for species',i))}
        }
      }
    }
  }
 }

