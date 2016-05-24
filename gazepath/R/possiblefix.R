possiblefix <-
function(speed, speed.cutoff){
  rle <- rle(ifelse(speed > speed.cutoff, 's', 'f'))
  return(rle$lengths[rle$values == 'f'])
}
