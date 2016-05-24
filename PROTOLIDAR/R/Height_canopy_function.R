Height_canopy_function <-
function(data_3D,distance_left,distance_right,min_canopy,max_canopy){

  x_plant <- y_plant <- z_plant <- NULL

  canopy <- subset(data_3D, data_3D$z_plant >= distance_left & data_3D$z_plant <= distance_right & data_3D$y_plant >= min_canopy & data_3D$y_plant <= max_canopy,select=c(x_plant,y_plant,z_plant))

  mean_height_canopy <- mean(canopy[,2])

  min_height_canopy  <- min(canopy[,2])

  max_height_canopy  <- max(canopy[,2])

  return(data.frame(mean_height_canopy,min_height_canopy,max_height_canopy))

}

