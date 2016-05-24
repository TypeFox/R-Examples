Sub.A <-
function(z, t, Mat, i, j){
if(t == 3){
index.a <- which(Mat[,i] == 1 & Mat[,j] == 0 & Mat[, -c(i, j)] == 0)
index.b <- which(Mat[,i] == 0 & Mat[,j] == 1 & Mat[, -c(i, j)] == 0)
index.ab <- which(Mat[,i] == 1 & Mat[,j] == 1 & Mat[, -c(i, j)] == 0)
Aij <- z[index.a] + z[index.b] + 2 * z[index.ab]
}

else {
index.a <- which(Mat[,i] == 1 & Mat[,j] == 0 & rowSums(Mat[, -c(i, j)]) == 0)
index.b <- which(Mat[,i] == 0 & Mat[,j] == 1 & rowSums(Mat[, -c(i, j)]) == 0)
index.ab <- which(Mat[,i] == 1 & Mat[,j] == 1 & rowSums(Mat[, -c(i, j)]) == 0)
Aij <- z[index.a] + z[index.b] + 2 * z[index.ab]
}
names(Aij) <- paste("A", i, j, sep="")
return(as.numeric(Aij))
}
