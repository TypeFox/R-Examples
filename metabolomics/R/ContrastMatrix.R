ContrastMatrix<-function(contrasts, levels){

  mat<-makeContrasts(contrasts=contrasts,levels=levels)
  return(mat)
}