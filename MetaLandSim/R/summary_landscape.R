summary_landscape <-
function(object)
  {
  
  if (class(object)!="landscape") 
  {
  stop(paste(object, " should be an object of class class 'landscape'.", sep=""), call. = FALSE)
  }
  
  land_area <- (object$mapsize^2)/10000
  n_patches <- object$number.patches
  mean_area <- object$mean.area
  SD_area <- object$SD.area
  
  m <- matrix.graph(object,"euc_distance")
  
  m2 <- upper.tri(m, diag = FALSE)
  
  m <- m*m2
  
  m[m == 0] <- NA
  
  mean_dist <- mean(m,na.rm=TRUE)
  min_dist <- min(m,na.rm=TRUE)
  
  output <- as.data.frame(matrix(ncol=1,nrow=6))
  
  rownames(output) <- c("landscape area (hectares)","number of patches","mean patch area (hectares)","SD patch area","mean distance amongst patches (meters)","minimum distance amongst patches (meters)")
  colnames(output) <- "Value"
  
  output[1,1] <- round(land_area,3)
  output[2,1] <- round(n_patches,3)
  output[3,1] <- round(mean_area,3)
  output[4,1] <- round(SD_area,3)
  output[5,1] <- round(mean_dist,3)
  output[6,1] <- round(min_dist,3)
  
  
  return(output)
  }
