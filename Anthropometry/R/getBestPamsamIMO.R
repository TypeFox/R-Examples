getBestPamsamIMO <- function(data,maxsplit,orness=0.7,type,ah,verbose,...){

 if (nrow(data) <= maxsplit){
  maxsplit <- max(nrow(data) - 1, 2)
 }
  
 if(maxsplit > 2){
  DIST <- ext.dist(data, maxsplit, orness, ah, verbose) 
  out <- INCAnumclu(DIST, K = maxsplit, method = "pam", L = NULL, noise = NULL)
   if(max(out$INCAindex[2:maxsplit]) <= 0.2){
    k <- 3
   }else{
     #According to Irigoien at al. (2008), "We propose to choose k as the value of k prior to the first biggest 
     #slope decrease".
     #out$INCAindex[2:maxsplit] are the out$INCAindex values without the NA value located at the first position.
     #The diff function calculates the subtraction between consecutive elements. Thus, we can know the decrease 
     #between INCA indexes for each k.
     restas <- diff(out$INCAindex[2:maxsplit], 1)  
     if(sum(restas >= 0) == length(restas)){#if all subtraction are positive, there are no decreases). 
     #In this radical case, we fix k = 3.
      k <- 3
     }else{
       #Greatest subtraction:
       min_resta <- min(restas)
       if(length(which(restas == min_resta)) > 1){ #it may happen that several subtractions take the same value. 
         #In this case, R will provide several k and will display many warnings on the screen. Therefore, we 
         #choose among those subtractions, the first one (because the others take the same value and therefore 
         #they are not smaller):
       k <- which(restas == min_resta)[1] + 1 #In order to identify correctly the value of k prior to the first 
       #biggest slope decrease, we have to add 1.
      }else{
        k <- which(restas == min_resta) + 1
       }
      }
     }
   x.ps <- pamsam(data, k = k, type = type, DIST = DIST, ah = ah, verbose = verbose, ...)
 }else{
   DIST <- ext.dist(data, maxsplit, orness, ah, verbose)
   x.ps <- pamsam(data, k = 2, type = type, DIST = DIST, ah = ah, verbose = verbose,...)
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
