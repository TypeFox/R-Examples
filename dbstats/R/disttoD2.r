

  ###########################
  #### function D2toDist ####
  ###########################


 disttoD2 <- function(distance){
   if (class(distance)[1]!="D2"){
     if (class(distance)[1]!="dist"&&class(distance)[1]!="dissimilarity") 
       stop("the distance matrix must be of class 'dist'/'dissimilarity")
     if (class(distance)[1]=="dissimilarity"&&attr(distance,"Metric")=="mixed")
      Delta <- as.matrix(distance)
     else
      Delta<-as.matrix(distance)^2
   }
   else
    Delta<-distance
   
   class(Delta)<-"D2"             
  return(Delta)
 }
