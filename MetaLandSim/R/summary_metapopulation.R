summary_metapopulation <-
function(object)
  {
  
  if (class(object)!="metapopulation") 
  {
  stop(paste(object, " should be an object of class class 'metapopulation'.", sep=""), call. = FALSE)
  }
  
  land_area <- (object$mapsize^2)/10000
  n_patches <- object$number.patches
  mean_area <- object$mean.area
  SD_area <- object$SD.area
  
  object1 <- object$nodes.characteristics
  
  n_snapshots <- ncol(object1)-8
  
  species_names <- paste(rep("species occurrence - snapshot",n_snapshots),as.character(1:n_snapshots))
  
  object2 <- remove.species(object)
  
  m <- matrix.graph(object2,"euc_distance")
  
  m2 <- upper.tri(m, diag = FALSE)
  
  m <- m*m2
  
  m[m == 0] <- NA
  
  mean_dist <- mean(m,na.rm=TRUE)
  min_dist <- min(m,na.rm=TRUE)
  
  output <- as.data.frame(matrix(ncol=1,nrow=6+n_snapshots))
  
  rownames(output) <- c("landscape area (hectares)","number of patches","mean patch area (hectares)","SD patch area","mean distance amongst patches (meters)","minimum distance amongst patches (meters)",species_names)
  colnames(output) <- "Value"
  
  output[1,1] <- land_area
  output[2,1] <- n_patches
  output[3,1] <- mean_area
  output[4,1] <- SD_area
  output[5,1] <- mean_dist
  output[6,1] <- min_dist
  
  for (i in 9:ncol(object1)){
		occ <- sum(object1[,i])
		occ_p <- (occ*100)/n_patches
		output[i-2,1] <- occ_p		
		}

  output[,1] <- round(output[,1],3)		
  
  return(output)
  
  }
