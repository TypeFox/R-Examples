Sub.S <-
function(z, t, Mat, j){
index <- which(Mat[,j] == 1 & rowSums(Mat[, -j]) == 0)
Sj <- sum(z[index])
names(Sj) <- paste("S", j, sep="")
return(Sj)
}
