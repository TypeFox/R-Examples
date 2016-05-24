biweight <- function(ty,M,c){
  if (ty<=c){
    w <- 0
  } else if (ty>c & ty<M){
    w <- (1-((ty-M)/(c-M))^2)^2
  } else if (ty>=M){
    w <- 1
  }
}