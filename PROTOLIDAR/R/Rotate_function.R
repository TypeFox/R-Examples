Rotate_function <-
function(data_3D,angle){

  z <- -data_3D[,3]

  x <-  data_3D[,1]

  x_rot <-c(x*cos(angle)-z*sin(angle))

  y_rot <-c(x*sin(angle)+z*cos(angle))

  return(data.frame(x_rot,y_rot))

}

