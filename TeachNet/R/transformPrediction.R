transformPrediction <- function(pred, threshold){
  if(pred<threshold){pred=0}else{pred=1}
  return(pred)
}