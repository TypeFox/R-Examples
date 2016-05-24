## ################################################################
## Function to find the exposure/treatment variable in a given data
## set
## ################################################################

find.treat <- function(data,
                       treat)
{
  if (missing(data)){
    stop("Argument 'object' contains no data.") 
  }else{  
      
    if(is.character(treat)){                  
      if(any(treat == names(data))){

        name.treat <- names(data[charmatch(treat, names(data))])
        treat      <- data[,charmatch(treat, names(data))]
        
      }else{
        stop("Argument 'treat' is not found in data.")
      }
      
    }else{

      if(is.numeric(treat)){
        if(length(treat) == 1){
          if(treat <= as.numeric(ncol(data))){
            
            name.treat   <- names(data)[treat]
            treat        <- data[,treat]
            
          }else stop("Argument 'treat' is not found in data.")
        }else stop("Argument 'treat' is not correctly given.")
      }
    }
  }

  if (length(levels(as.factor(treat))) != 2)
    stop("Argument 'treat' has either only one value or more than two values.")
  
  return(list(treat      = treat,
              name.treat = name.treat))
  
}



 




