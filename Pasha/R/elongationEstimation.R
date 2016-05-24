

elongationEstimation <- function(piledValuesPlus, piledValuesMinus, stepShift) 
{
    #void elongationEstimation(int *piledValuesPlus, int *piledValuesMinus, int *stepShift, int *maxShift, int *piledValues_Size, int *stepShift_Size) 
    
    # Padding with 0s if one vector is longer than the other one
    maxLength <- max(c(length(piledValuesPlus),length(piledValuesMinus)))
    
    length(piledValuesPlus) <- maxLength
    length(piledValuesMinus) <- maxLength
    
    piledValuesPlus[is.na(piledValuesPlus)] <- 0
    piledValuesMinus[is.na(piledValuesMinus)] <- 0
    
    functionReturn <- .C("C_elongationEstimation", as.integer(piledValuesPlus), as.integer(piledValuesMinus), res = as.double(stepShift), as.integer(max(stepShift)), as.integer(maxLength), as.integer(length(stepShift)))
    
    return(functionReturn$res)
}

#elongationEstimationOld <- function(piledValuesPlus, piledValuesMinus, stepShift) 
#{
#    # Padding with 0s if one vector is longer than the other one
#    maxLength=max(c(length(piledValuesPlus),length(piledValuesMinus)))
#
#    length(piledValuesPlus)=maxLength
#    length(piledValuesMinus)=maxLength
#
#    piledValuesPlus[is.na(piledValuesPlus)]=0
#    piledValuesMinus[is.na(piledValuesMinus)]=0
#
#    # Computing the overlap score for each shifting step
#    res=sapply(stepShift, function(currentShift)
#    {
#                return(sum(c(rep(0,currentShift),piledValuesPlus)*c(piledValuesMinus,rep(0,currentShift))))
#                #return(sum(piledValuesPlus[1:(length(piledValuesPlus-currentShift))]*piledValuesMinus[currentShift:length(piledValuesMinus)]))
#    })
#
#    return(res)
#}
