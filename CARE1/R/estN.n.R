estN.n <-
function(z){
t <- as.integer(log(length(z) + 1) / log(2))
Mat <- matrix(0, ncol = t, nrow = 2 ^ t - 1)
for(i in 1:(2 ^ t - 1)) Mat[i,] <- rev(as.integer(intToBits(i))[1:t])
size <- sapply(1:t, function(i) Sub.n(z, t, Mat, i))
return(size)
}
