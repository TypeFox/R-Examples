postvar <-
  function(Kmat,Jmat,Cisqrt){
  ### get posterior variance
  
  pvar <- rowSums((cbind(Kmat,Jmat)%*%Cisqrt)^2)
  
}
