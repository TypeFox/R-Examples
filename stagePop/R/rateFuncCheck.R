#' rateFunCheck
#' 
#' check the user-defined rate functions exist and are in the correct format
#'
#' checks existence of deathFunc, reproFunc, durationFunc, immigrationFunc and develFunc if 
#' there is time dependent duration.
#'
#' Also checks that output from rate functions is a single value
#'
#' @param x Vector of state variables within the DDE solver. To access a variable use: x$speciesName['stageName',strainNumber]
#' @param numSpecies Number of species
#' @param numStages Number of life stages
#' @param numStrains Number of strains for each species
#' @param timeDependDuration vector of logicals defining whether the stage duration is time dependent for each species
#' @param rateFunctions List of rate functions to check
#' @param speciesNames Vector of species names
#' @param stageNames Vector of species' stage names
#' @param startTime First entry in timeVec vector
#;
#' @return TRUE/FALSE - TRUE if rate functions are OK, FALSE if not.
rateFuncCheck=function(x,numSpecies,numStages,numStrains,timeDependDuration,rateFunctions,speciesNames,stageNames,startTime){

#  print('enter rateFuncCheck')
  
  flag=TRUE
  
  if (! is.function(rateFunctions$deathFunc)) {
    warning('Please define deathFunc(stage,x,time,species,strain)')
    flag=FALSE
  }

  if (! is.function(rateFunctions$reproFunc)) {
    warning('Please define reproFunc(x,time,species,strain)')
    flag=FALSE
  }

  if (! is.function(rateFunctions$durationFunc)) {
    warning('Please define durationFunc(stage,x,time,species,strain)')
    flag=FALSE
  }

  if (! is.function(rateFunctions$immigrationFunc)) {
    warning('Please define immigrationFunc(stage,x,time,species)')
    flag=FALSE
  }

  if (! is.function(rateFunctions$emigrationFunc)) {
    warning('Please define emigrationFunc(stage,x,time,species,strain)')
    flag=FALSE
  }
  
  if ( sum(timeDependDuration) > 0 ){
    if (! is.function(rateFunctions$develFunc)) {
      warning('Please define develFunc(stage,x,time,species,strain)')
      flag=FALSE
    }
  }


  if (flag){
    
    for (i in seq(1,numSpecies)){
      
      df=rateFunctions$reproFunc(x,startTime,i,1)

      if (is.nan(df)){print(paste('Output from reproFunc is not defined for species',i)); flag=FALSE}else if (length(df)>1){print('Output from reproFunc must be a single value (not a vector)'); flag=FALSE}
      
      for (j in seq(1,max((numStages[i]-1),1))){
        
        df=rateFunctions$deathFunc(j,x,startTime,i,1)
        if (is.nan(df)){print(paste('Output from deathFunc is not defined for species',i,'stage',j)); flag=FALSE}else if (length(df)>1){print('Output from deathFunc must be a single value (not a vector)'); flag=FALSE}
        
      
        if (numStages[i]>1){
          df=rateFunctions$durationFunc(j,x,startTime,i,1)
          if (is.nan(df)){print(paste('Output from durationFunc is not defined for species',i,'stage',j)); flag=FALSE}else if (length(df)>1){print('Output from durationFunc must be a single value (not a vector)'); flag=FALSE}  
        }
        
        if (timeDependDuration[i]){
          df=rateFunctions$develFunc(j,x,startTime,i,1)
          if (is.nan(df)){print(paste('Output from develFunc is not defined for species',i,'stage',j)); flag=FALSE}else if (length(df)>1){print('Output from develFunc must be a single value (not a vector)'); flag=FALSE}    
        }
      }
    }
  }
#  print('leave rateFuncCheck')

  return(flag)

}
