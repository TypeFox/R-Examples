get.incremental<-function(triangle)
{
  triangle<-as.matrix(triangle)  
  triangle.inc<-triangle
  m<-nrow(triangle)
  for (i in 1:m) triangle.inc[i,1:(m-i+1)]<-c(as.numeric(triangle[i,1]),diff(as.numeric(triangle[i,1:(m-i+1)])))
  return(triangle.inc)
}