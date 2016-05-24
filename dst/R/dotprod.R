dotprod<-function(SR,BD,g,f){
  # Generalized inner product of two matrix
  # 
  if (ncol(SR) !=nrow(BD)) {
    stop("nb cols mat1 not equal nb rows mat2")
  }
  ff<-match.fun(f); gg<-match.fun(g)
  resul<-matrix(FALSE,nrow(SR),ncol(BD))
  for(i in 1:nrow(SR)){
    for(j in 1:ncol(BD)){
      temp<-ff(SR[i,1],BD[1,j])
        for (k in 2:ncol(SR)) {
          temp<-gg(temp,ff(SR[i,k],BD[k,j]))
        }
      resul[i,j]=temp
    }
  }
  resul
}
