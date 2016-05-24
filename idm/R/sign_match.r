sign_match<-function(matA,matB){
  
  p=ncol(matA)
  na=nrow(matA)
  nb=nrow(matB)
  n=min(na,nb)
  A.sign <- sign(matA[1:n,])
  B.sign <- sign(matB[1:n,])
  s.Mat <- matrix(rep(1,(n*p)),n,p)
  s.Mat[A.sign != B.sign] = -1
  
  s.tr=sign(apply(s.Mat,2,sum))
  #this fixes the case when s.tr colsum equals to zero
  s.tr[which(s.tr==0)]=-1
  matB = matB %*% diag(s.tr)
  #   Bg.sign=  sign(matB[1,])
  #   print("Bg")
  #   print(Bg.sign)
  
  matB
  
}