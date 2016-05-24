
PDCDistance <- function(x, y, ...){
  tryCatch({
  as.numeric(pdcDist(cbind(x, y), ...))},
  error=function(e){print(e); NA})   
}

