PCanalysis <-function(TFBS, nPCs, Sequences){

  model<-pca(TFBS, nPcs=nPCs, method="svd", center=TRUE)
  matrix.residuals<-residuals(model,Sequences, nPcs=model@nPcs)
  residus<-apply(matrix.residuals,1,function(vector){vector%*%vector})
}

