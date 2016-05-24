pvcam <-
function(xba, batch, y, threshold=0.6) {

  pvcaresult <- mypvcaBatchAssess(xba, factordata=data.frame(batch=batch, y=y), threshold=threshold)
  pvcavec <- as.vector(pvcaresult$dat)
  names(pvcavec) <- pvcaresult$label
  pvcalist <- as.list(pvcavec)
  names(pvcalist)[1] <- "batch_y"
  # 'pvcalist' is a list with the following element names: "batch_y", "batch", "y", "resid"
  
  return(pvcalist$y)

}
