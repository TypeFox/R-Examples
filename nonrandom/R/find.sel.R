
## ######################################################
## Funtion to find variables selected in a given data set
## ######################################################

find.sel <- function(data,
                     sel,
                     sel.name="sel")
{
  
  if(missing(data) | is.null(data)){
    stop("Argument 'object' contains no data.")      
  }else{    
    if(is.data.frame(sel)){      
      seldata <- sel      
    }else{
      if(is.character(sel)){
        if ( is.na(sum(charmatch(sel, names(data)))) ){ 
          stop(paste("One or more entries in argument '",
                     sel.name,
                     "' are not contained in data.", sep=""))
        }else{
          if ( any(charmatch(sel, names(data)) == 0)  ){             
            stop(paste("One or more entries in argument '",
                       sel.name,
                       "' are in line with multiple variables in data.", sep=""))
          }else{
            seldata <- as.data.frame(data[,charmatch(sel, names(data))])
            names(seldata) <- sel          
          }
        } 
      }else{        
        if(is.numeric(sel)){
          if (max(sel) <= ncol(data)){            
            seldata <- as.data.frame(data[, sel])
            names(seldata) <- names(data)[sel]               
          }else{
            stop(paste("One or more entries in argument '",
                       sel.name,
                       "' are not contained in data.", sep=""))  
          }       
        }
      }
    }
  }
  return(seldata)  
}




