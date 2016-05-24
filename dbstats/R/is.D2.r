  ########################
  #### function is.D2 ####
  ########################


 is.D2 <- function(x){
   aux<-ifelse(class(x)=="D2",TRUE,FALSE)
   return(aux)
 }