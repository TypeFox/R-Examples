getBestPamsamMO <- function(data,maxsplit,orness = 0.7,type,ah,verbose, ...){
 
 if(nrow(data) <= maxsplit){
  maxsplit <- max(nrow(data) - 1, 2)
 }
  
 if(maxsplit > 2){
  #Create the pamsams with 2 to maxsplit splits:
  max.asw <- -1
  for (k in 2 : maxsplit){
   xtemp.ps <- pamsam(data, k = k, type = type, DIST = NULL, maxsplit = maxsplit, 
                      orness = orness, ah, verbose, ...)
   if (xtemp.ps$asw > max.asw){
    max.asw <- xtemp.ps$asw
    x.ps <- xtemp.ps
   }
  }   
 }else{
   x.ps <- pamsam(data, k = 2, type = type, DIST = NULL, maxsplit = maxsplit, 
                  orness = orness, ah, verbose, ...)
  }
 
 #Object checking:
 check <- c()
 for(i in 1:length(x.ps)){
  check[i] <- exists(names(x.ps)[i], where = x.ps)  
 }
  
 if("FALSE" %in% check){
  stop("Any of the objects doesn't exist!. Revise the function.") 
 }else{
  return(x.ps) 
 }
}
