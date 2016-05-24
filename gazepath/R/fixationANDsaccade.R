fixationANDsaccade <-
function(speed, thres_vel, thres_dur, Hz){
  thres_dur <- thres_dur * (Hz / 1000)
  fixsac <- ifelse(speed > thres_vel, 's', 'f')
  ## Set a minimum saccade duration of 10 ms
  rle <- rle(fixsac)
  fixsac[cumsum(rle$lengths)[which(rle$lengths < 10 * (Hz/1000) & rle$values == 's')]] <- 'f'
  
  while(length(which(rle$lengths < 10 * (Hz/1000) & rle$values == 's')) != 0){
    rle <- rle(fixsac)
    fixsac[cumsum(rle$lengths)[which(rle$lengths < 10 * (Hz/1000) & rle$values == 's')]] <- 'f'
  }
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
