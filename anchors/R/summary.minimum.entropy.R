summary.minimum.entropy <- function(object, ..., average=FALSE) {
  ## calculate the proportion in each \hat{C}
  if (average) 
    return(tapply( object$Prop, object$MinEnt, sum))
  else
    return( as.data.frame(object))
}
