"chop" <-
function(y,z,off){
#  Chops off leading and trailing NAs from a vector

  n <- length(y)
  if( !is.na(y[1])& !is.na(y[n])){
    list(y=y,z=z,off=off)}
  else{
    if( is.na(y[1])) Recall(y[-1],z[-1,],off[-1])
    else Recall(y[-n],z[-n,],off[-n])
  }
}

