stdz <- function(x, weight=NULL){
  if(is.null(weight)){
    weight <- rep(1, length(x))
  }
  x <- x-wtd.mean(x, weight, na.rm=TRUE)
  x <- x/sqrt(wtd.var(x, weight, na.rm=TRUE))
  x
}
