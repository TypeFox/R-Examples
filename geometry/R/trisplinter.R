"trisplinter" <-
function(T,p,threshold=sqrt(.Machine$double.eps)){
  rownorm2 = function(x) drop(sqrt((x^2)%*%c(1,1)))
  d1 = p[T[,1],] - p[T[,2],]
  d2 = p[T[,2],] - p[T[,3],]
  d1 = d1 / rownorm2(d1)
  d2 = d2 / rownorm2(d2)
  ar = d1[,1]*d2[,2] - d1[,2]*d2[,1]
  return(abs(ar) < threshold)
}
