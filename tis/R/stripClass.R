stripClass <- function(x, classString){
  classVec <- class(x)
  newClassVec <- classVec[classVec != classString]
  if(length(newClassVec) == 0) return(unclass(x))
  else {
    class(x) <- newClassVec
    return(x)
  }
}
