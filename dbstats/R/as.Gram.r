  ##########################
  #### function as.Gram ####
  ##########################


 as.Gram <- function(x){
   G<-as.matrix(x)
   class(G)<-"Gram"
   return(G)
 }