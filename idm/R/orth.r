orth<- function(A){
  ## questa function produce una base ortogonale Q per il range di A.
  ## Dunque t(Q) %*% Q =I, quindi le colonne di Q si trovano nello stesso
  ## spazio delle colonne di A. Il numero di colonne di Q equivale al rango di A.
  
  n=nrow(A)
  p=ncol(A)
  
  Asvd=svd(A)
  r=qr(A)$rank
  Q=Asvd$u[,1:r]
  Q
}