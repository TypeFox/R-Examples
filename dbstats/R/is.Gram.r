  ##########################
  #### function is.Gram ####
  ##########################


 is.Gram <- function(x){
   aux<-ifelse(class(x)=="Gram",TRUE,FALSE)
   return(aux)
 }