Extract_plant_3D_function <-
function(out,z_min,z_max,y_min,y_max,distance_left,distance_right){ 

  x_cm <- y_cm <- z_cm <- NULL

  data_stem <- subset(out,out$y_cm > y_min & out$y_cm<= y_max & out$z_cm >= z_min  & out$z_cm <= z_max,select=c(x_cm,y_cm,z_cm))

  x_c       <- out$x_cm - min(data_stem$x_cm)

  y_c       <- out$y_cm

  z_c       <- out$z_cm - min(data_stem$z_cm)

  data_cero <- data.frame(x_c,y_c,z_c)

  data_plant <- subset(data_cero, data_cero$z_c >= distance_left & data_cero$z_c <= distance_right,select=c(x_c,y_c,z_c))

  x_plant <- data_plant[,1]

  y_plant <- data_plant[,2]

  z_plant <- data_plant[,3]

  return(data.frame(x_plant,y_plant,z_plant))

  }

