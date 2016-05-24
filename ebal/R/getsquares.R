getsquares <-
function(mat){
mat.add <- matrix(NA,nrow(mat),0)
mat.add.names <- c()
for(i in 1:ncol(mat)){
 if(dim(table(mat[,i]))>2)
  {
 mat.add <- cbind(mat.add,mat[,i]^2)
 mat.add.names <- c(paste(colnames(mat)[i],".2",sep=""),mat.add.names)
 }
}
colnames(mat.add) <- mat.add.names
out <- cbind(mat,mat.add)
}

