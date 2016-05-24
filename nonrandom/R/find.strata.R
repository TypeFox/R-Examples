find.strata <- function(data,     ## a data frame
                        strata,   ## a string or numeric
                        match)  ## TRUE if match is given
{

  if (missing(data)){

    stop("Argument 'object' contains no data.") 
    
  }else{  
      
    if(is.character(strata)){                  
      if(any(strata == names(data))){

        name.stratum.index   <- names(data[charmatch(strata, names(data))])
        stratum.index        <- data[,charmatch(strata, names(data))]
        levels.stratum.index <- levels(as.factor(stratum.index))
        
      }else{
        stop("Argument 'stratum.index' is not found in data.")
      }      
    }else{

      if(is.numeric(strata)){      
        if(length(strata) == 1){
          if(strata <= ncol(data)){

            name.stratum.index   <- names(data)[strata]
            stratum.index        <- data[,strata]
            levels.stratum.index <- levels(as.factor(stratum.index))
            
          }else stop("Argument 'stratum.index' is not found in data.")
        }else{
          stop("Argument 'stratum.index' contains two or more entries.")
        }
      }
    }
  }

  if (!match)
    if(any(levels(as.factor(stratum.index)) == "0"))
      stop("Argument 'stratum.index' may not involve '0' as stratum index.")

  return(list(stratum.index        = stratum.index,
              name.stratum.index   = name.stratum.index,
              levels.stratum.index = levels.stratum.index))
  
}
