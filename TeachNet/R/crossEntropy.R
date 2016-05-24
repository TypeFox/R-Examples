crossEntropy <- function(x,y) {
  if(y==1) {return(-log(x))}
  else {return(-log(1-x))}
}