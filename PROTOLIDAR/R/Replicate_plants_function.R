Replicate_plants_function <-
function(plants_rotate,data_3D,latitude,longitude){

  x_rot <- plants_rotate[,1]

  y_rot <- plants_rotate[,2]

  z     <- data_3D[,2]

  rep_z <- rep(z,length(latitude)) 

  rep_X <- rep(x_rot,length(latitude))

  rep_Y <- rep(y_rot,length(latitude))



  dup_xcoord <- rep(latitude ,each=length(x_rot))

  dup_ycoord <- rep(longitude,each=length(y_rot))



  XCOORD <- rep_X + dup_xcoord

  YCOORD <- rep_Y + dup_ycoord

  return (data.frame(XCOORD,YCOORD,z))

  }

