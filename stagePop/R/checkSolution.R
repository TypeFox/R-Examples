#' checkSolution
#' 
#' Check whether any of the state variables are significantly negative and emit suitable warnings.
#'
#' @param output Model output from \code{\link{popModel}}
#' @param numSpecies Number of species
#' @param numStages Number of life stages
#' @param numStrains Number of strains for each species
#' @param ntol Negative tolerance value (i.e. a warning is produced if variable<-(ntol*max(variable))). The larger ntol is, the larger the negative values that are tolerated. 
#'
#' @seealso \code{\link{popModel}}
#'
#' @return Nothing if there are no problems with the output, otherwise warnings are generated
#' @export
checkSolution=function(output,numSpecies,numStages,numStrains,ntol){
  flag=FALSE
  
  ct=2
  for (i in seq(1,numSpecies)){
    for (k in seq(1,numStrains[i])){
      for (j in seq(1,numStages[i])){
        data=output[,ct]
        if (length(data[data<(-ntol*max(data))])>0){
          warning(paste('There are negative values for',names(output[1,ct]),'The minimum value is',min(data)))
          flag=TRUE}
        ct=ct+1
      }
    }
  }
  
  if (flag) {
    warning('Decreasing the tolerance values in the DDE solver may help remove negative values. If this does not work then check for errors in your rate functions and other popModel inputs.')
  }
}

        
