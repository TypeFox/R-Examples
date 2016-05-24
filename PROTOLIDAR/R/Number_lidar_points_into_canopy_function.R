Number_lidar_points_into_canopy_function <-
function (data_3D,distance_left,distance_right,min_canopy,max_canopy){

  x_plant <- y_plant <- z_plant <- NULL

  canopy <- subset(data_3D, data_3D$z_plant >= distance_left & data_3D$z_plant <= distance_right & data_3D$y_plant >= min_canopy & data_3D$y_plant <= max_canopy,select=c(x_plant,y_plant,z_plant))

  N_points <- length(canopy[,1])

  return (N_points)

  }

