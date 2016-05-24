arrMat<-function(input){
  d<-dim(input)
  if(length(d)==3){mat <- matrix(input, nrow = 1, byrow = TRUE)}
  if(length(d)==4){mat <- matrix(input, nrow = d[4], byrow = TRUE)}
  return(mat)
}