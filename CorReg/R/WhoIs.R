# ' Give the partition implied by a structure
# ' 
# ' @param Z the structure (square matrix)
# ' @param I1,I2,I3 wanted outputs (I1= explicative covariates, I2=redundant covariates, I3=independent covariates)
# ' 
# '@export 
WhoIs<-function(Z=Z,I3=F,I2=T,I1=T){
  res=list()
  if(I2){
    res$I2=which(colSums(Z)!=0)# Qui est ? gauche
  }
  if(I1){
    res$I1=which(colSums(Z)==0)#qui n'est pas ? gauche
  }
  if(I3){
    res$I3=which(colSums(Z)==0 & rowSums(Z)==0)# qui est totalement isol?
  }
  return(res)
}