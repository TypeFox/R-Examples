which.distance.function <- function(method){
  if(method=='levenshtein'){
    distance.function<-levenshtein.distance
  }
  if(method=='hamming'){
    distance.function<-hamming.distance
  }
  if(method=='levenshtein.damerau'){
    distance.function<-levenshtein.damerau.distance
  }
  return(distance.function)
}