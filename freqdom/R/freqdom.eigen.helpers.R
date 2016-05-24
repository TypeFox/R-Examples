freqdom.inflambdas = function(E){
  apply(E$values,2,min)
}
freqdom.infalphas = function(E){
  D = dim(E$values)
  
  R = abs(E$values[,1:(D[2]-1)] - E$values[,1 + 1:(D[2]-1)])
  
  A = R[,1:(D[2]-2)]
  B = R[,1 + 1:(D[2]-2)]
  R = array(0,c(dim(A),2))
  R[,,1] = A
  R[,,2] = B
  
  R = apply(R,1:2,min)
  R = cbind(A[,1],R,B[,D[2]-2])
  
  apply(R,2,min)
}
freqdom.invlamintegral = function(E){
  V = 1/E$values
  1/apply(V,2,mean)
}
