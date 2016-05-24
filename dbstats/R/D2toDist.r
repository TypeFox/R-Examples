

  ###########################
  #### function D2toDist ####
  ###########################


 D2toDist <- function(D2){

  if (class(D2)!="D2")
    stop(" 'D2' must be of class D2")
    
  distance <-as.dist(D2^.5)
  return(distance)
 }
 