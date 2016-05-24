generateu <-
function(mat){
  names=colnames(mat)
  r=length(names)
  U=list(r)
  for(i in 1:r){
    U[[i]]=IndMatrix(mat[,i])
    colnames(U[[i]])=paste(names[i],"(",colnames(U[[i]]),")",sep="")
  }
  names(U)=names
  class(U)="Umatrix"
  return(U)
}
