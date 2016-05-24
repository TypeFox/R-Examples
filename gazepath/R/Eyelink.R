Eyelink <-
function(X, Y, D, height_mm, width_mm, height_px, width_px, Hz, defliction = .1, thres_dur = 100){
  speed <- Speed_Deg(X, Y, D, height_mm, width_mm, height_px, width_px, Hz)
  accel <- speed^2
  distance <- c(.001, sqrt(diff(X)^2 + diff(Y)^2))
  #transform .1 deg to number of pixels
  offset <- tan((defliction / 2) * pi/180) * D * (1 / (height_mm / height_px)) * 2
  
  output <- ifelse(speed > 30 & accel > 8000 & distance > offset, 's', 'f')
  output[is.na(speed)] <- NA
  if(Hz > 250){
    for(i in 1:(Hz / 250)){
      index <- which(rle(output)$values == 's' & rle(output)$lengths <= 2)
      output[cumsum(rle(output)$lengths)[index]] <- 'f'
    }
  }
  ## use a duration threshold
  thres_dur <- thres_dur * (Hz / 1000)
  rle <- rle(output)
  classify <- numeric()
  for(i in 1:length(rle$values)){
    if(is.na(rle$values[i])){
      classify <- c(classify, rep(NA, rle$lengths[i]))
    } else{
      if(rle$values[i] == 'f' & rle$lengths[i] >= thres_dur){
        classify <- c(classify, rep('f', rle$lengths[i]))
      }
      if(rle$values[i] == 'f' & rle$lengths[i] < thres_dur){
        classify <- c(classify, rep('u', rle$lengths[i]))
      }
      if(rle$values[i] == 's'){
        classify <- c(classify, rep('s', rle$lengths[i]))
      }
    }
  }
  return(classify)
}
