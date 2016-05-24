LWA_lidar_function <-
function(height_canopy,ground_area,row_spacing){

 LWA <- 2* height_canopy * (ground_area/row_spacing)

return(LWA)

}

