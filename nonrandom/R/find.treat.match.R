## ################################################################
## Function to find the exposure/treatment variable in a given data
## set for matching; difference to find.treat: only one (but
## different) value of treat is allowed per object
## ################################################################

find.treat.match <- function(data,
                             treat,
                             obj.name)
{
  if ( missing(data) ){
    stop(paste("Argument 'object.",obj.name,"' contains no data.",
               sep=""))
  }else{  
      
    if( is.character(treat) ){                  
      if( any(treat == names(data)) ){

        name.treat <- names(data[charmatch(treat, names(data))])
        treat      <- data[,charmatch(treat, names(data))]
        
      }else{
        stop(paste("Argument 'treat' is not found in data of 'object.",obj.name,".",
                   sep=""))
      }
      
    }else{

      if( is.numeric(treat) ){
        if( length(treat) == 1 ){
          if( treat <= as.numeric(ncol(data)) ){
            
            name.treat   <- names(data)[treat]
            treat        <- data[,treat]
            
          }else stop(paste("Argument 'treat' is not found in data in data of 'object.",
                           obj.name,".", sep=""))
        }else stop(paste("Argument 'treat' is not correctly given in data of 'object.",
                         obj.name,".", sep=""))
      }
    }
  }

  if ( length(levels(as.factor(treat))) != 1 )
    stop(paste("Argument 'treat' has more than one value in object.",
               obj.name,".", sep=""))
  
  return(list(treat      = treat,
              name.treat = name.treat))
  
}



 




