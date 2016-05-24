## This function calculates the smooth for the precision function
MFW <- function(x){
  if(sum(is.na(x)) / length(x) > .5){
    return(rep(NA, length(x)))
  } else {
    return(rep(median(x, na.rm = T), length(x)))
  }
}