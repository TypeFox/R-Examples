estN.pair <-
function(z){
t <- as.integer(log(length(z) + 1) / log(2))
Mat <- matrix(0, ncol = t, nrow = 2 ^ t - 1)
for(i in 1:(2 ^ t - 1)) Mat[i,] <- rev(as.integer(intToBits(i))[1:t])
estN.ij <- Sub.pair(z, t, Mat, 1, 2)
for(i in 1:(t - 1)){
for(j in (i + 1):t){
estN.ij = rbind(estN.ij, Sub.pair(z, t, Mat, i, j))
}
}
return(estN.ij[-1,])
}
