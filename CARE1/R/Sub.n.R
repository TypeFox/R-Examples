Sub.n <-
function(z, t, Mat, j){
index <- which(Mat[,j] == 1)
nj <- sum(z[index])
names(nj) <- paste("n", j, sep="")
return(nj)
}
