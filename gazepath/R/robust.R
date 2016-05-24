robust <-
function(x, Hz) {
  mis <- rle(ifelse(is.na(x), 0, 1))
  return(mean(mis$length[mis$values == 1]) / Hz)
}
