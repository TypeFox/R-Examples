TRV_lidar_function <-
function(height_canopy,width_canopy,row_spacing){

TRV <-height_canopy * width_canopy * 10000  / row_spacing

return(TRV)

}

