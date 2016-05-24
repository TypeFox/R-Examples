
# Dynamic time warping distance is calculated by using dtw package

DTWDistance <- function(x, y, ...){
  tryCatch({
  dtw(x, y, ...)$distance}, 
  error=function(e){print(e); NA})   
}