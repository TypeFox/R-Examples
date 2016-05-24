Sub.B <-
function(z, t, Mat, i, j){
index <- which(Mat[,i] == 1 & Mat[,j] == 1)
Bij <- sum(z[index])
names(Bij) <- paste("B", i, j, sep="")
return(Bij)
}
