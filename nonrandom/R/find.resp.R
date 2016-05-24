## ##########################################################
## Function to find the response variable in a give ndata set
## ##########################################################

find.resp <- function(data,
                      resp)
{
  
  if (missing(data)){
    stop("Argument 'object' contains no data.") 
  }else{  
      
    if(is.character(resp)){

      if(any(resp == names(data))){

        name.resp <- names(data[charmatch(resp, names(data))])
        resp      <- data[,charmatch(resp, names(data))]
        
      }else{
        stop("Argument 'resp' is not found in data.")
      }      
    }else{

      if(is.numeric(resp)){
        if(length(resp) == 1){
          if(resp <= as.numeric(ncol(data))){
            
            name.resp <- names(data)[resp]
            resp      <- data[,resp]
            
          }else stop("Argument 'resp' is not found in data.")
        }else stop("Argument 'resp' is not correctly given.")
      }
    }
  }
  
  return(list(resp      = resp,
              name.resp = name.resp))
  
}



 




