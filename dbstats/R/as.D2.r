

  ########################
  #### function as.D2 ####
  ########################


 as.D2 <- function(x){
 
  # if (class(x)=="D2")
#     Delta <- x
#   else{  
#    if (class(x)[1]=="dissimilarity"&&attr(x,"Metric")=="mixed")
#      Delta <- as.matrix(x)
#    else
#      Delta <- as.matrix(x)^2
#   }
#

   Delta = as.matrix(x)
   class(Delta)<-"D2"
 
   return(Delta)
 }
