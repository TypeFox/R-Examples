# function standardizes everything to 0 - 1
range01 <- function(x){
  if(min(x) == max(x)){
    x <- rep(0, length(x))
    return(x)
  }
  (x-min(x))/(max(x)-min(x))
}
