FrechetDistance <- function (x, y, tx=NULL, ty=NULL, ...) {
  
  if(is.null(tx)){
    tx <- c(1:length(x))
  }
  
  if(is.null(ty)){
    ty <- c(1:length(y))
  }
  
  tryCatch({
  distFrechet(tx, x, ty, y, ...)},
  error=function(e) {print(e); NA})   
}



