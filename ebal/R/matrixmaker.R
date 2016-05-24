matrixmaker <-
function(mat){
k <- ncol(mat)
obs <- nrow(mat)
out <- matrix(NA,obs,((k*(k+1))/2 +1))
count <- 0
 for (i in 1:k){
  for (j in 1:i){
  count <- count + 1
   out[,count] <- mat[,i] * mat[,j]
   colnames(out)[c(count,((k*(k+1))/2 +1))] <-
   c(paste(colnames(mat)[i],".",colnames(mat)[j],sep=""), "dummy")
  }
 }
out <- cbind(mat,out[,1:((k*(k+1))/2)])
out
}

